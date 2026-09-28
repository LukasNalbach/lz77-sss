/**
 * part of LukasNalbach/lz77-sss
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include <lz77_sss/misc/utils.hpp>

template <typename value_t, typename pred_t>
inline uint64_t parallel_partition(
    std::vector<value_t>& values, uint64_t beg, uint64_t end, uint16_t p, pred_t pred)
{
    using segment_t = std::pair<uint64_t, uint64_t>;
    const uint64_t size = end - beg;
    auto chunk_beg = [&](uint64_t i_p) { return beg + size * i_p / p; };
    std::vector<uint64_t> split(p);

    #pragma omp parallel for num_threads(p)
    for (uint16_t i_p = 0; i_p < p; i_p++) {
        split[i_p] = std::partition(
            values.begin() + chunk_beg(i_p),
            values.begin() + chunk_beg(i_p + 1),
            pred) - values.begin();
    }

    uint64_t mid = beg;

    for (uint16_t i_p = 0; i_p < p; i_p++) {
        mid += split[i_p] - chunk_beg(i_p);
    }

    std::vector<segment_t> left_segs;
    std::vector<segment_t> right_segs;
    uint64_t num_misplaced = 0;

    for (uint16_t i_p = 0; i_p < p; i_p++) {
        const uint64_t left_end = std::min(chunk_beg(i_p + 1), mid);
        const uint64_t right_beg = std::max(chunk_beg(i_p), mid);

        if (split[i_p] < left_end) {
            left_segs.emplace_back(split[i_p], left_end);
            num_misplaced += left_end - split[i_p];
        }

        if (right_beg < split[i_p]) {
            right_segs.emplace_back(right_beg, split[i_p]);
        }
    }

    if (num_misplaced == 0) return mid;

    auto seek = [](const std::vector<segment_t>& segs, uint64_t offset) {
        uint64_t seg = 0;

        while (offset >= segs[seg].second - segs[seg].first) {
            offset -= segs[seg].second - segs[seg].first;
            seg++;
        }

        return std::make_pair(seg, segs[seg].first + offset);
    };

    #pragma omp parallel for num_threads(p)
    for (uint16_t i_p = 0; i_p < p; i_p++) {
        uint64_t from = num_misplaced * i_p / p;
        const uint64_t to = num_misplaced * (i_p + 1) / p;
        if (from == to) continue;

        auto [left_seg, left_pos] = seek(left_segs, from);
        auto [right_seg, right_pos] = seek(right_segs, from);

        while (from < to) {
            const uint64_t len = std::min({ to - from,
                left_segs[left_seg].second - left_pos,
                right_segs[right_seg].second - right_pos });

            std::swap_ranges(
                values.begin() + left_pos,
                values.begin() + left_pos + len,
                values.begin() + right_pos);

            from += len;
            left_pos += len;
            right_pos += len;

            if (left_pos == left_segs[left_seg].second && ++left_seg < left_segs.size()) {
                left_pos = left_segs[left_seg].first;
            }

            if (right_pos == right_segs[right_seg].second && ++right_seg < right_segs.size()) {
                right_pos = right_segs[right_seg].first;
            }
        }
    }

    return mid;
}

template <typename value_t, typename key_fnc_t>
inline void parallel_nth_element(
    std::vector<value_t>& values, uint64_t beg, uint64_t nth, uint64_t end, uint16_t p, key_fnc_t key)
{
    constexpr uint64_t min_par_chunk_size = 1 << 14;
    constexpr uint64_t max_sample_size = 1 << 14;
    std::mt19937_64 gen;
    std::vector<uint64_t> sample;

    while (true) {
        const uint64_t size = end - beg;
        const uint16_t p_used = std::min<uint64_t>(p, size / min_par_chunk_size);
        if (p_used <= 1) break;

        const uint64_t sample_size = std::min<uint64_t>(max_sample_size, std::sqrt(size));
        const double quantile = (nth - beg) / double(size);
        const uint64_t target = quantile * sample_size;
        const uint64_t offset = 1 + 3 * std::sqrt(sample_size * quantile * (1 - quantile));
        const bool keep_left = 2 * (nth - beg) < size;
        const uint64_t pivot_idx = keep_left ?
            std::min(target + offset, sample_size - 1) :
            target - std::min(target, offset);

        std::uniform_int_distribution<uint64_t> pos_distrib(beg, end - 1);
        sample.resize(sample_size);

        for (uint64_t& value : sample) {
            value = key(values[pos_distrib(gen)]);
        }

        std::nth_element(sample.begin(), sample.begin() + pivot_idx, sample.end());
        const uint64_t pivot = sample[pivot_idx];
        auto less = [&](const value_t& value) { return key(value) < pivot; };
        auto less_equal = [&](const value_t& value) { return key(value) <= pivot; };

        if (keep_left) {
            const uint64_t mid = parallel_partition(values, beg, end, p_used, less);

            if (nth < mid) {
                end = mid;
            } else {
                beg = parallel_partition(values, mid, end, p_used, less_equal);
                if (nth < beg) return;
            }
        } else {
            const uint64_t mid = parallel_partition(values, beg, end, p_used, less_equal);

            if (nth >= mid) {
                beg = mid;
            } else {
                end = parallel_partition(values, beg, mid, p_used, less);
                if (nth >= end) return;
            }
        }
    }

    std::nth_element(values.begin() + beg, values.begin() + nth, values.begin() + end,
        [&](const value_t& a, const value_t& b) { return key(a) < key(b); });
}
