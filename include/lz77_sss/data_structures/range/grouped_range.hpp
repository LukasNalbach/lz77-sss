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
#include <cstdint>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

#include <lz77_sss/data_structures/bit_aligned_interleaved_vectors.hpp>
#include <lz77_sss/data_structures/range/range.hpp>
#include <lz77_sss/misc/utils.hpp>

template <typename impl_t>
class grouped_range final : public range_ds {
    static constexpr uint64_t min_par_chunk_size = 1 << 15;

    uint64_t num_points = 0;
    std::vector<uint64_t> G_S;
    std::vector<impl_t> R_g;

    inline uint64_t num_groups() const { return R_g.size(); }

    inline uint64_t frq(uint64_t g) const { return G_S[g + 1] - G_S[g]; }

    inline uint64_t group_of(uint64_t coord) const
    {
        return uint64_t(std::upper_bound(G_S.begin() + 1, G_S.end() - 1, coord) - G_S.begin()) - 1;
    }

    inline void to_internal(uint64_t g, point_t& point) const
    {
        point.x = uint64_t(point.x) - G_S[g];
        point.y = uint64_t(point.y) - G_S[g];
    }

    inline void to_external(uint64_t g, point_t& point) const
    {
        point.x = uint64_t(point.x) + G_S[g];
        point.y = uint64_t(point.y) + G_S[g];
    }

    inline point_t internal_point(const std::vector<uint32_t>& grp,
        const bit_aligned_interleaved_vectors<3>& P, uint64_t i) const
    {
        point_t point { .x = P.get<0>(i), .y = P.get<1>(i), .weight = P.get<2>(i) };
        to_internal(grp[i], point);
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

    void build_groups(const std::vector<uint32_t>& chr, const bit_aligned_interleaved_vectors<3>& P, uint16_t p)
    {
        const uint64_t n = chr.size();
        std::vector<uint32_t> chr_by_y;
        no_init_resize(chr_by_y, n);

        #pragma omp parallel for num_threads(p) schedule(static)
        for (uint64_t i = 0; i < n; i++) chr_by_y[P.get<1>(i)] = chr[i];

        std::vector<std::vector<uint64_t>> starts(p);

        #pragma omp parallel num_threads(p)
        {
            const uint64_t i_p = omp_get_thread_num();
            const uint64_t q = omp_get_num_threads();

            for (uint64_t y = std::max<uint64_t>(1, n * i_p / q); y < n * (i_p + 1) / q; y++) {
                if (chr_by_y[y] != chr_by_y[y - 1]) starts[i_p].emplace_back(y);
            }
        }

        chr_by_y = std::vector<uint32_t>();
        const uint64_t limit = std::max<uint64_t>(1, n / max_group_share_inv);
        G_S.assign(1, 0);
        uint64_t prev = 0;

        auto add_block = [&](uint64_t end) {
            const uint64_t beg = prev;
            prev = end;
            if (beg == G_S.back()) return;
            if (end - G_S.back() > limit) G_S.emplace_back(beg);
        };

        for (const std::vector<uint64_t>& local : starts) {
            for (uint64_t y : local) add_block(y);
        }

        add_block(n);
        G_S.emplace_back(n);
    }

    void build_buckets(const std::vector<uint32_t>& grp,
        const bit_aligned_interleaved_vectors<3>& P, uint16_t p)
    {
        const uint64_t n = grp.size();
        const uint64_t m = num_groups();
        const uint16_t q = std::max<uint64_t>(1, std::min<uint64_t>(p, n / min_par_chunk_size));
        std::vector<uint64_t> offs(q * m, 0);

        #pragma omp parallel num_threads(q)
        {
            const uint64_t i_p = omp_get_thread_num();
            uint64_t* cnt = offs.data() + i_p * m;
            for (uint64_t i = n * i_p / q; i < n * (i_p + 1) / q; i++) cnt[grp[i]]++;
        }

        for (uint64_t g = 0; g < m; g++) {
            uint64_t pos = G_S[g];

            for (uint16_t i_p = 0; i_p < q; i_p++) {
                const uint64_t cnt = offs[i_p * m + g];
                offs[i_p * m + g] = pos;
                pos += cnt;
            }
        }

        bit_aligned_vector idx(n, n);

        #pragma omp parallel num_threads(q)
        {
            const uint64_t i_p = omp_get_thread_num();
            uint64_t* nxt = offs.data() + i_p * m;
            for (uint64_t i = n * i_p / q; i < n * (i_p + 1) / q; i++) idx.set_parallel(nxt[grp[i]]++, i);
        }

        auto build_bucket = [&](uint64_t g, uint16_t p_g) {
            std::vector<point_t> bucket;
            no_init_resize(bucket, frq(g));

            #pragma omp parallel for num_threads(p_g) schedule(dynamic, 65536)
            for (uint64_t k = 0; k < frq(g); k++) {
                bucket[k] = internal_point(grp, P, idx[G_S[g] + k]);
            }

            R_g[g] = impl_t(bucket, frq(g), p_g);
        };

        std::vector<uint64_t> small_buckets;

        for (uint64_t g = 0; g < m; g++) {
            if (frq(g) == 0) continue;

            if (frq(g) >= std::max<uint64_t>(min_par_chunk_size, n / p)) {
                build_bucket(g, p);
            } else {
                small_buckets.emplace_back(g);
            }
        }

        #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
        for (uint64_t k = 0; k < small_buckets.size(); k++) {
            build_bucket(small_buckets[k], 1);
        }
    }

public:
    static constexpr uint64_t max_group_share_inv = 256;

    grouped_range(const std::vector<uint32_t>& chr, const bit_aligned_interleaved_vectors<3>& P, uint16_t p)
    {
        const uint64_t n = chr.size();
        build_groups(chr, P, p);
        R_g.resize(G_S.size() - 1);
        std::vector<uint32_t> grp;
        no_init_resize(grp, n);

        #pragma omp parallel for num_threads(p) schedule(static)
        for (uint64_t i = 0; i < n; i++) grp[i] = uint32_t(group_of(P.get<1>(i)));

        if constexpr (impl_t::is_static()) num_points = n;

        if constexpr (std::is_same_v<impl_t, static_weighted_square_grid>) {
            if (sorted_by_weight(P, p)) {
                const uint64_t win_size = bench_win_size != 0 ? bench_win_size : impl_t::default_win_size;
                std::vector<uint64_t> extents(num_groups());
                for (uint64_t g = 0; g < num_groups(); g++) extents[g] = frq(g);

                impl_t::build(std::span(R_g), std::span<const uint64_t>(extents), n, win_size, p,
                    [&](uint64_t i) { return uint64_t(grp[i]); },
                    [&](uint64_t i) { return internal_point(grp, P, i); });

                return;
            }
        }

        build_buckets(grp, P, p);
    }

    bool is_decomposed() const override { return true; }
    bool is_static() const override { return impl_t::is_static(); }
    bool is_dynamic() const override { return impl_t::is_dynamic(); }
    std::string name() const override { return "d-" + impl_t::name(); }
    uint64_t size() const override { return num_points; }

    uint64_t size_in_bytes() const override
    {
        uint64_t bytes = sizeof(*this) + G_S.size() * sizeof(uint64_t);
        for (uint64_t g = 0; g < num_groups(); g++) bytes += R_g[g].size_in_bytes();
        return bytes;
    }

    void insert(uint64_t, point_t point) override
    {
        if constexpr (impl_t::is_dynamic()) {
            const uint64_t g = group_of(point.y);
            to_internal(g, point);
            R_g[g].insert(point);
            num_points++;
        }
    }

    result_t lighter_point_in_range(uint64_t, uint64_t weight,
        uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const override
    {
        if constexpr (impl_t::is_static()) {
            if (x1 > x2 || y1 > y2) return { point_t { }, false };
            const uint64_t g = group_of(y1);
            const uint64_t beg = G_S[g];
            auto [point, found] = R_g[g].lighter_point_in_range(weight, x1 - beg, x2 - beg, y1 - beg, y2 - beg);
            to_external(g, point);
            return { point, found };
        } else {
            return { point_t { }, false };
        }
    }

    result_t point_in_range(uint64_t,
        uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const override
    {
        if constexpr (impl_t::is_dynamic()) {
            if (x1 > x2 || y1 > y2) return { point_t { }, false };
            const uint64_t g = group_of(y1);
            const uint64_t beg = G_S[g];
            auto [point, found] = R_g[g].point_in_range(x1 - beg, x2 - beg, y1 - beg, y2 - beg);
            to_external(g, point);
            return { point, found };
        } else {
            return { point_t { }, false };
        }
    }
};
