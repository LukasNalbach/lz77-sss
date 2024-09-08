#pragma once

#include <lz77_sss/data_structures/sample_index/sample_index.hpp>

template <typename pos_t, typename sidx_t, typename lce_r_t>
template <direction dir>
void sample_index<pos_t, sidx_t, lce_r_t>::build_samples(pos_t typ_lce_r, bool log) {
    auto time = now();

    for (uint32_t c = 0; c < (1 << 16); c++) {
        SXIV2[dir][c].b = no_occ;
        SXIV2[dir][c].e = no_occ;
    }

    if (log) {
        std::cout << "building SLC" << (dir == LEFT ? "S" : "P") << std::flush;
    }

    std::vector<pos_t> SLCX;
    SLCX.reserve(s + 1);
    SLCX.emplace_back(0);

    for (sidx_t i = 1; i < s; i++) {
        pos_t pos_im1 = S[SXA<dir>(i - 1)];
        pos_t pos_i = S[SXA<dir>(i)];
        pos_t lcx = lce<dir>(pos_im1, pos_i, max_lce_l);
        SLCX.emplace_back(lcx);
    }

    SLCX.emplace_back(0);

    if (log) {
        time = log_runtime(time);
        std::cout << "sorting SLC" << (dir == LEFT ? "S" : "P") << std::flush;
    }

    uint64_t baseline = malloc_count_current();
    std::vector<pos_t> SLCX_sorted = SLCX;
    ips4o::sort(SLCX_sorted.begin(), SLCX_sorted.end() - 1);

    if (log) {
        time = log_runtime(time);
    }

    sidx_t slcx_rng_min = bin_search_min_geq<pos_t, sidx_t>(
        3, 0, s - 1, [&](sidx_t i){return SLCX_sorted[i];});

    if (typ_lce_r == std::numeric_limits<pos_t>::max()) {
        typ_lce_r = SLCX_sorted[s - 1];
    }

    pos_t max_sampled_len = std::min<pos_t>(
        SLCX_sorted[s - 1],
        dir == LEFT ? max_lce_l : typ_lce_r);

    sidx_t slcx_rng_max = bin_search_min_geq<pos_t, sidx_t>(
        max_sampled_len, 0, s - 1, [&](sidx_t i){return SLCX_sorted[i];});

    double max_num_samples = 2.0 * s;
    double slcx_rng = slcx_rng_max - slcx_rng_min;
    smpl_pat_lens[dir] = {1, 2};
    SXIVX<dir>().resize(2);
    int16_t num_pat_lens = 2 + std::floor(
        (max_num_samples - slcx_rng / 2.0) /
        (double{slcx_rng_min} + slcx_rng / 2.0)
    );

    for (int16_t i = 0; i < num_pat_lens; i++) {
        double rel_slcx_rnk = (i + 1.0) / double{num_pat_lens};
        pos_t slcx_rnk = slcx_rng_min + rel_slcx_rnk * slcx_rng;
        pos_t len = SLCX_sorted[slcx_rnk];
        while (len <= smpl_pat_lens[dir].back()) len++;
        if (len > max_sampled_len) break;

        smpl_pat_lens[dir].emplace_back(len);
        slcx_rnk = bin_search_min_geq<pos_t, sidx_t>(
            len, 0, s - 1, [&](sidx_t i){
                return SLCX_sorted[i];
        });

        SXIVX<dir>().emplace_back(sxivx_map_t<dir>(0,
            sxivx_hash<dir>(this, len),
            sxivx_eq<dir>(this, len)
        ));

        SXIVX<dir>().back().reserve(slcx_rnk);
    }

    num_pat_lens = smpl_pat_lens[dir].size();
    std::vector<sidx_t> sxiv_b(num_pat_lens, 0);
    SLCX_sorted.clear();
    SLCX_sorted.shrink_to_fit();

    if (log) {
        std::cout << "chose " << num_pat_lens - 2 << " pattern"
            << " lengths in the range [" << smpl_pat_lens[dir][2]
            << ", " << smpl_pat_lens[dir].back() << "]";
        time = log_runtime(time);
        std::cout << "sampling S" << (dir == LEFT ? "P" : "S")
            << "A intervals" << std::flush;
    }

    for (sidx_t i = 1; i <= s; i++) {
        pos_t lcx = SLCX[i];

        for (int16_t idx = num_pat_lens - 1; idx >= 2; idx--) {
            pos_t len = smpl_pat_lens[dir][idx];
            if (lcx >= len) break;
            pos_t pos_im1 = S[SXA<dir>(i - 1)];

            if (is_pos_in_T<dir>(pos_im1, len - 1)) {
                SXIVX<dir>()[idx].emplace(
                    sxa_interval_t {sxiv_b[idx], i - 1});
            }

            sxiv_b[idx] = i;
        }

        if (lcx < 2) {
            pos_t pos_im1 = S[SXA<dir>(i - 1)];

            if (is_pos_in_T<dir>(pos_im1, 1)) {
                uint16_t pat = val_offs<uint16_t, dir, 0>(T, pos_im1);
                SXIV2[dir][pat].b = sxiv_b[1];
                SXIV2[dir][pat].e = i - 1;
            }

            sxiv_b[1] = i;

            if (dir == LEFT && lcx < 1) {
                uint8_t pat = val_offs<uint8_t, RIGHT, 0>(T, pos_im1);
                SCIV[pat].b = sxiv_b[0];
                SCIV[pat].e = i - 1;
                sxiv_b[0] = i;
            }
        }
    }

    #ifndef NDEBUG
    for (uint8_t idx = 0; idx < num_pat_lens; idx++) {
        sxiv_b[idx] = 0;
    }

    for (sidx_t i = 1; i <= s; i++) {
        pos_t lcx = SLCX[i];

        for (int16_t idx = num_pat_lens - 1; idx >= 2; idx--) {
            pos_t len = smpl_pat_lens[dir][idx];
            if (lcx >= len) break;
            pos_t pos_im1 = S[SXA<dir>(i - 1)];

            if (is_pos_in_T<dir>(pos_im1, len - 1)) {
                pos_hash = pos_im1;
                auto it = SXIVX<dir>()[idx].find({s, s});
                assert(it != SXIVX<dir>()[idx].end());
                assert(it->b == sxiv_b[idx]);
                assert(it->e == i - 1);
            }
            
            sxiv_b[idx] = i;
        }

        if (lcx < 2) {
            pos_t pos_im1 = S[SXA<dir>(i - 1)];

            if (is_pos_in_T<dir>(pos_im1, 1)) {
                uint16_t pat = val_offs<uint16_t, dir, 0>(T, pos_im1);
                assert(SXIV2[dir][pat].b == sxiv_b[1]);
                assert(SXIV2[dir][pat].e == i - 1);
            }

            sxiv_b[1] = i;

            if (dir == LEFT && lcx < 1) {
                uint8_t pat = val_offs<uint8_t, RIGHT, 0>(T, pos_im1);
                assert(SCIV[pat].b == sxiv_b[0]);
                assert(SCIV[pat].e == i - 1);
                sxiv_b[0] = i;
            }
        }
    }
    #endif

    if (log) {
        std::cout << " (" << format_size(malloc_count_current() - baseline) << ")";
        time = log_runtime(time);
    }
}