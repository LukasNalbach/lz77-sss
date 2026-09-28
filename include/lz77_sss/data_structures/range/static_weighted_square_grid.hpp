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
#include <span>
#include <vector>

#include <lz77_sss/data_structures/bit_aligned_interleaved_vectors.hpp>
#include <lz77_sss/data_structures/range/range_point.hpp>
#include <lz77_sss/misc/utils.hpp>

class static_weighted_square_grid {
public:
    static constexpr bool is_decomposed() { return false; }
    static constexpr bool is_static() { return true; }
    static constexpr bool is_dynamic() { return false; }

    static constexpr uint64_t default_win_size = 16384;

    using point_t = range_point_t;

protected:
    struct __attribute__((packed)) window {
        uint40_t beg;
        uint16_t len;
    };

    static constexpr uint64_t min_par_chunk_size = 1 << 15;

    uint64_t num_points = 0;
    uint64_t win_size = 0;
    uint64_t grid_width = 0;
    uint64_t pos_max_excl = 0;
    bit_aligned_interleaved_vectors<3> points;
    std::vector<window> windows;

    inline uint64_t grid_index(uint64_t x_w, uint64_t y_w) const { return grid_width * y_w + x_w; }
    inline uint64_t window_index(point_t point) const { return grid_index(point.x / win_size, point.y / win_size); }

    inline point_t point_at(uint64_t i) const
    {
        return point_t { .x = points.get<0>(i), .y = points.get<1>(i), .weight = points.get<2>(i) };
    }

public:
    static_weighted_square_grid() = default;

    static_weighted_square_grid(
        std::vector<point_t>& input_points,
        uint64_t pos_max_excl, uint16_t p = 1,
        uint64_t win_size = default_win_size)
    {
        if (bench_win_size != 0) win_size = bench_win_size;
        bool sorted = true;

        #pragma omp parallel for num_threads(p) reduction(&& : sorted)
        for (uint64_t i = 1; i < input_points.size(); i++) {
            sorted = sorted && input_points[i - 1].weight <= input_points[i].weight;
        }

        if (!sorted) {
            ips4o::parallel::sort(input_points.begin(), input_points.end(),
                [](const point_t& p1, const point_t& p2) { return p1.weight < p2.weight; }, p);
        }

        build(std::span(this, 1), std::array<uint64_t, 1> { pos_max_excl }, input_points.size(), win_size, p,
            [](uint64_t) { return uint64_t { 0 }; },
            [&](uint64_t i) { return input_points[i]; });

        input_points.clear();
        input_points.shrink_to_fit();
    }

    template <typename grid_of_t, typename point_at_t>
    static void build(
        std::span<static_weighted_square_grid> grids, std::span<const uint64_t> extents,
        uint64_t num_points, uint64_t win_size, uint16_t p,
        grid_of_t grid_of, point_at_t input_point)
    {
        const uint64_t num_grids = grids.size();
        const uint64_t max_weight = num_points == 0 ? 0 : uint64_t(input_point(num_points - 1).weight);
        std::vector<uint64_t> win_offs(num_grids + 1, 0);

        for (uint64_t g = 0; g < num_grids; g++) {
            grids[g].win_size = win_size;
            grids[g].pos_max_excl = extents[g];
            grids[g].grid_width = div_ceil(extents[g], win_size);
            win_offs[g + 1] = win_offs[g] + grids[g].grid_width * grids[g].grid_width;
        }

        const uint64_t num_win = win_offs[num_grids];
        const uint16_t q = std::max<uint64_t>(1, std::min<uint64_t>({ p,
            num_points / min_par_chunk_size, num_points / (2 * num_win + 1) }));
        std::vector<uint64_t> offs(q * num_win, 0);

        #pragma omp parallel num_threads(q)
        {
            const uint64_t i_p = omp_get_thread_num();
            uint64_t* cnt = offs.data() + i_p * num_win;

            for (uint64_t i = num_points * i_p / q; i < num_points * (i_p + 1) / q; i++) {
                const uint64_t g = grid_of(i);
                cnt[win_offs[g] + grids[g].window_index(input_point(i))]++;
            }
        }

        std::vector<uint64_t> word_offs(num_grids + 1, 0);

        for (uint64_t g = 0; g < num_grids; g++) {
            static_weighted_square_grid& grid = grids[g];
            const uint64_t num_win_g = win_offs[g + 1] - win_offs[g];
            no_init_resize(grid.windows, num_win_g);
            uint64_t pos = 0;

            for (uint64_t w = 0; w < num_win_g; w++) {
                const uint64_t beg = pos;

                for (uint64_t i_p = 0; i_p < q; i_p++) {
                    uint64_t& off = offs[i_p * num_win + win_offs[g] + w];
                    const uint64_t cnt = off;
                    off = pos;
                    pos += cnt;
                }

                grid.windows[w] = window { .beg = beg, .len = uint16_t(pos - beg) };
            }

            grid.num_points = pos;

            if (pos != 0) {
                grid.points = bit_aligned_interleaved_vectors<3>(
                    pos, { grid.pos_max_excl, grid.pos_max_excl, max_weight }, init_mode::none);
            }

            word_offs[g + 1] = word_offs[g] + grid.points.num_words();
        }

        #pragma omp parallel num_threads(q)
        {
            const uint64_t i_p = omp_get_thread_num();
            const uint64_t words_beg = word_offs[num_grids] * i_p / q;
            const uint64_t words_end = word_offs[num_grids] * (i_p + 1) / q;
            uint64_t g = std::upper_bound(word_offs.begin(), word_offs.end(), words_beg) - word_offs.begin() - 1;

            for (; g < num_grids && word_offs[g] < words_end; g++) {
                const uint64_t from = std::max(words_beg, word_offs[g]);
                const uint64_t to = std::min(words_end, word_offs[g + 1]);
                if (from < to) grids[g].points.zero_words(from - word_offs[g], to - word_offs[g]);
            }

            #pragma omp barrier

            uint64_t* nxt = offs.data() + i_p * num_win;

            for (uint64_t i = num_points * i_p / q; i < num_points * (i_p + 1) / q; i++) {
                const uint64_t g = grid_of(i);
                static_weighted_square_grid& grid = grids[g];
                const point_t point = input_point(i);
                const uint64_t pos = nxt[win_offs[g] + grid.window_index(point)]++;
                grid.points.init_parallel(pos, { point.x, point.y, point.weight });
            }
        }
    }

    inline uint64_t size() const { return num_points; }

    inline uint64_t size_in_bytes() const
    {
        return sizeof(*this) + windows.size() * sizeof(window) + points.size_in_bytes();
    }

    std::tuple<point_t, bool> lighter_point_in_range(
        uint64_t weight,
        uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const
    {
        uint64_t xw_1 = x1 / win_size;
        uint64_t xw_2 = x2 / win_size;
        uint64_t yw_1 = y1 / win_size;
        uint64_t yw_2 = y2 / win_size;
 
        uint64_t xi_beg = xw_1 + (x1 % win_size != 0);
        uint64_t yi_beg = yw_1 + (y1 % win_size != 0);
        uint64_t xi_end = xw_2 + (x2 % win_size == (win_size - 1));
        uint64_t yi_end = yw_2 + (y2 % win_size == (win_size - 1));
        bool has_inner_windows = xi_beg < xi_end && yi_beg < yi_end;
 
        if (has_inner_windows) {
            uint64_t w_idx = grid_index(xi_beg, yi_beg);
            uint64_t nxt_row_offs = grid_width - (xi_end - xi_beg);
 
            for (uint64_t y_w = yi_beg; y_w < yi_end; y_w++) {
                for (uint64_t x_w = xi_beg; x_w < xi_end; x_w++) {
                    const window& w = windows[w_idx];

                    if (w.len != 0) {
                        const point_t point = point_at(w.beg);
                        if (point.weight < weight) return { point, true };
                    }
                        
                    w_idx++;
                }
 
                w_idx += nxt_row_offs;
            }
        }

        uint64_t w_idx = grid_index(xw_1, yw_1);
        uint64_t nxt_row_offs = grid_width - (xw_2 - xw_1 + 1);
 
        for (uint64_t y_w = yw_1; y_w <= yw_2; y_w++) {
            bool y_contained = has_inner_windows && yi_beg <= y_w && y_w < yi_end;
            uint64_t x_w = xw_1;
 
            while (x_w <= xw_2) {
                if (y_contained && x_w == xi_beg) {
                    w_idx += xi_end - xi_beg;
                    x_w = xi_end;
                    continue;
                }
 
                const window& w = windows[w_idx];
                uint64_t end = w.beg + w.len;
 
                for (uint64_t i = w.beg; i < end; i++) {
                    const point_t point = point_at(i);
                    if (point.weight >= weight) break;
 
                    if (x1 <= point.x && point.x <= x2 &&
                        y1 <= point.y && point.y <= y2
                    ) return { point, true };
                }
 
                w_idx++;
                x_w++;
            }
 
            w_idx += nxt_row_offs;
        }
 
        return { { 0, 0, 0 }, false };
    }

    static constexpr std::string name() { return "swsg"; }
};