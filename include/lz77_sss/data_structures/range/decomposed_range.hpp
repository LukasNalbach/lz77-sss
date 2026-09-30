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

#include <array>
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

#include <lz77_sss/data_structures/bit_aligned_interleaved_vectors.hpp>
#include <lz77_sss/data_structures/range/range.hpp>
#include <lz77_sss/misc/utils.hpp>

template <typename impl_t>
class decomposed_range final : public range_ds {
    static constexpr uint64_t min_par_chunk_size = 1 << 15;

    using offsets_t = std::vector<std::array<uint64_t, 256>>;

    uint64_t num_points = 0;
    std::array<uint64_t, 257> C_S = { 0 };
    std::array<impl_t, 256> R_c;

    inline uint64_t frq(uint8_t c) const { return C_S[c + 1] - C_S[c]; }

    inline void to_internal(uint8_t c, point_t& point) const
    {
        uint64_t rnk_c = C_S[c];
        point.x = uint64_t(point.x) - rnk_c;
        point.y = uint64_t(point.y) - rnk_c;
    }

    inline void to_internal(uint8_t c,
        uint64_t& x1, uint64_t& x2, uint64_t& y1, uint64_t& y2) const
    {
        uint64_t rnk_c = C_S[c];
        x1 -= rnk_c;
        x2 -= rnk_c;
        y1 -= rnk_c;
        y2 -= rnk_c;
    }

    inline void to_external(uint8_t c, point_t& point) const
    {
        uint64_t rnk_c = C_S[c];
        point.x = uint64_t(point.x) + rnk_c;
        point.y = uint64_t(point.y) + rnk_c;
    }

    inline point_t internal_point(const std::vector<uint8_t>& chr,
        const bit_aligned_interleaved_vectors<3>& P, uint64_t i) const
    {
        point_t point { .x = P.get<0>(i), .y = P.get<1>(i), .weight = P.get<2>(i) };
        to_internal(chr[i], point);
        return point;
    }

    static bool sorted_by_weight(const bit_aligned_interleaved_vectors<3>& P, uint16_t p)
    {
        bool sorted = true;

        #pragma omp parallel for num_threads(p) reduction(&& : sorted)
        for (uint64_t i = 1; i < P.size(); i++) {
            sorted = sorted && P.get<2>(i - 1) <= P.get<2>(i);
        }

        return sorted;
    }

    void build_buckets(const std::vector<uint8_t>& chr,
        const bit_aligned_interleaved_vectors<3>& P, uint16_t p, offsets_t& offs)
    {
        const uint64_t n = chr.size();
        const uint16_t q = offs.size();
        bit_aligned_vector idx(n, n);

        #pragma omp parallel num_threads(q)
        {
            const uint64_t i_p = omp_get_thread_num();
            std::array<uint64_t, 256>& nxt = offs[i_p];

            for (uint64_t i = n * i_p / q; i < n * (i_p + 1) / q; i++) {
                idx.set_parallel(nxt[chr[i]]++, i);
            }
        }

        auto build_bucket = [&](uint8_t c, uint16_t p_c) {
            std::vector<point_t> bucket;
            no_init_resize(bucket, frq(c));

            #pragma omp parallel for num_threads(p_c) schedule(dynamic, 65536)
            for (uint64_t k = 0; k < frq(c); k++) {
                bucket[k] = internal_point(chr, P, idx[C_S[c] + k]);
            }

            R_c[c] = impl_t(bucket, frq(c), p_c);
        };

        std::vector<uint8_t> small_buckets;

        for (uint16_t c = 0; c < 256; c++) {
            if (frq(c) == 0) continue;

            if (frq(c) >= std::max<uint64_t>(min_par_chunk_size, n / p)) {
                build_bucket(c, p);
            } else {
                small_buckets.emplace_back(c);
            }
        }

        #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
        for (uint64_t k = 0; k < small_buckets.size(); k++) {
            build_bucket(small_buckets[k], 1);
        }
    }

public:
    decomposed_range(const std::vector<uint8_t>& chr,
        const bit_aligned_interleaved_vectors<3>& P, uint16_t p)
    {
        const uint64_t n = chr.size();
        const uint16_t q = std::max<uint64_t>(1, std::min<uint64_t>(p, n / min_par_chunk_size));
        offsets_t offs(q);

        #pragma omp parallel num_threads(q)
        {
            const uint64_t i_p = omp_get_thread_num();
            std::array<uint64_t, 256>& cnt = offs[i_p];
            cnt.fill(0);

            for (uint64_t i = n * i_p / q; i < n * (i_p + 1) / q; i++) {
                cnt[chr[i]]++;
            }
        }

        for (uint16_t c = 0; c < 256; c++) {
            uint64_t pos = C_S[c];

            for (uint16_t i_p = 0; i_p < q; i_p++) {
                const uint64_t cnt = offs[i_p][c];
                offs[i_p][c] = pos;
                pos += cnt;
            }

            C_S[c + 1] = pos;
        }

        if constexpr (impl_t::is_static()) num_points = n;

        if constexpr (std::is_same_v<impl_t, static_weighted_square_grid>) {
            if (sorted_by_weight(P, p)) {
                const uint64_t win_size = bench_win_size != 0 ? bench_win_size : impl_t::default_win_size;
                std::array<uint64_t, 256> extents;
                for (uint16_t c = 0; c < 256; c++) extents[c] = frq(c);

                impl_t::build(R_c, extents, n, win_size, p,
                    [&](uint64_t i) { return chr[i]; },
                    [&](uint64_t i) { return internal_point(chr, P, i); });

                return;
            }
        }

        build_buckets(chr, P, p, offs);
    }

    bool is_decomposed() const override { return true; }
    bool is_static() const override { return impl_t::is_static(); }
    bool is_dynamic() const override { return impl_t::is_dynamic(); }
    std::string name() const override { return "d-" + impl_t::name(); }
    uint64_t size() const override { return num_points; }

    uint64_t size_in_bytes() const override
    {
        uint64_t bytes = sizeof(*this);
        for (uint16_t c = 0; c < 256; c++) if (frq(c) != 0) bytes += R_c[c].size_in_bytes();
        return bytes;
    }

    void insert(uint64_t chr, point_t point) override
    {
        if constexpr (impl_t::is_dynamic()) {
            const uint8_t c = uint8_t(chr);
            to_internal(c, point);
            R_c[c].insert(point);
            num_points++;
        }
    }

    result_t lighter_point_in_range(uint64_t chr, uint64_t weight,
        uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const override
    {
        if constexpr (impl_t::is_static()) {
            const uint8_t c = uint8_t(chr);
            to_internal(c, x1, x2, y1, y2);
            auto [point, found] = R_c[c].lighter_point_in_range(weight, x1, x2, y1, y2);
            to_external(c, point);
            return { point, found };
        } else {
            return { point_t { }, false };
        }
    }

    result_t point_in_range(uint64_t chr,
        uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const override
    {
        if constexpr (impl_t::is_dynamic()) {
            const uint8_t c = uint8_t(chr);
            to_internal(c, x1, x2, y1, y2);
            auto [point, found] = R_c[c].point_in_range(x1, x2, y1, y2);
            to_external(c, point);
            return { point, found };
        } else {
            return { point_t { }, false };
        }
    }
};
