#pragma once

#include <lz77_sss/data_structures/sample_index/sample_index.hpp>

template <typename pos_t, typename sidx_t, typename lce_r_t>
void sample_index<pos_t, sidx_t, lce_r_t>::build_sxa12_intervals(uint16_t p, bool log) {
    auto time = now();

    if (log) {
        std::cout << "precomputing SPA1/2 and SSA1/2 intervals" << std::flush;
    }

    #pragma omp parallel for num_threads(p)
    for (uint16_t c = 0; c < 256; c++) {
        SCIV[c].b = bin_search_min_geq<bool, sidx_t>(
            true, 0, s, [&](sidx_t i){
                return i == s || char_to_uchar(T[S[SSA[i]]]) >= c;
        });

        if (SCIV[c].b == s || char_to_uchar(T[S[SSA[SCIV[c].b]]]) != c) {
            SCIV[c].b = no_occ;
            SCIV[c].e = no_occ;
        } else {
            SCIV[c].e = bin_search_max_lt<bool, sidx_t>(
                true, SCIV[c].b, s - 1, [&](sidx_t i){
                    return char_to_uchar(T[S[SSA[i]]]) > c;
            });
        }
    }

    for (uint16_t c1 = 0; c1 < 256; c1++) {
        #pragma omp parallel for num_threads(p)
        for (uint16_t c2 = 0; c2 < 256; c2++) {
            uint16_t c_l = (c1 << 8) | c2;
            uint16_t c_r = (c2 << 8) | c1;

            if (SCIV[c1].b == no_occ) {
                SXIV2[LEFT][c_l].b = no_occ;
                SXIV2[LEFT][c_l].e = no_occ;
                SXIV2[RIGHT][c_r].b = no_occ;
                SXIV2[RIGHT][c_r].e = no_occ;
            } else {
                SXIV2[LEFT][c_l].b = bin_search_min_geq<bool, sidx_t>(
                    true, SCIV[c1].b, SCIV[c1].e + 1, [&](sidx_t i){
                        return i == SCIV[c1].e + 1 || (S[SPA[i]] != 0 && char_to_uchar(T[S[SPA[i]] - 1]) >= c2);
                });

                if (SXIV2[LEFT][c_l].b == SCIV[c1].e + 1 || char_to_uchar(T[S[SPA[SXIV2[LEFT][c_l].b]] - 1]) != c2) {
                    SXIV2[LEFT][c_l].b = no_occ;
                    SXIV2[LEFT][c_l].e = no_occ;
                } else {
                    SXIV2[LEFT][c_l].e = bin_search_max_lt<bool, sidx_t>(
                        true, SXIV2[LEFT][c_l].b, SCIV[c1].e, [&](sidx_t i){
                            return S[SPA[i]] != 0 && char_to_uchar(T[S[SPA[i]] - 1]) > c2;
                    });
                }

                SXIV2[RIGHT][c_r].b = bin_search_min_geq<bool, sidx_t>(
                    true, SCIV[c1].b, SCIV[c1].e + 1, [&](sidx_t i){
                        return i == SCIV[c1].e + 1 || (S[SSA[i]] != n - 1 && char_to_uchar(T[S[SSA[i]] + 1]) >= c2);
                });

                if (SXIV2[RIGHT][c_r].b == SCIV[c1].e + 1 || char_to_uchar(T[S[SSA[SXIV2[RIGHT][c_r].b]] + 1]) != c2) {
                    SXIV2[RIGHT][c_r].b = no_occ;
                    SXIV2[RIGHT][c_r].e = no_occ;
                } else {
                    SXIV2[RIGHT][c_r].e = bin_search_max_lt<bool, sidx_t>(
                        true, SXIV2[RIGHT][c_r].b, SCIV[c1].e, [&](sidx_t i){
                            return S[SSA[i]] != n - 1 && char_to_uchar(T[S[SSA[i]] + 1]) > c2;
                    });
                }
            }
        }
    }

    if (log) {
        time = log_runtime(time);
    }
}

template <typename pos_t, typename sidx_t, typename lce_r_t>
template <direction dir>
void sample_index<pos_t, sidx_t, lce_r_t>::build_samples(pos_t typ_lce_r, uint16_t p, bool log) {
    auto time = now();

    if (log) {
        std::cout << "building SLC" << (dir == LEFT ? "S" : "P") << std::flush;
    }

    std::vector<pos_t> SLCX;
    no_init_resize(SLCX, s + 1);
    SLCX[0] = 0;
    SLCX[s] = 0;

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 1; i < s; i++) {
        SLCX[i] = lce<dir>(
            S[SXA<dir>(i - 1)],
            S[SXA<dir>(i)],
            max_lce_l
        );
    }

    if (log) {
        time = log_runtime(time);
        std::cout << "sorting SLC" << (dir == LEFT ? "S" : "P") << std::flush;
    }

    uint64_t baseline = malloc_count_current();
    std::vector<pos_t> SLCX_sorted = SLCX;
    ips4o::parallel::sort(SLCX_sorted.begin(), SLCX_sorted.end() - 1);

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

    std::vector<sidx_t> pat_len_rank(num_pat_lens, 0);

    for (int16_t i = 2; i < num_pat_lens; i++) {
        double rel_slcx_rnk = (double{i} + 1.0) / double{num_pat_lens};
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

        pat_len_rank[i] = slcx_rnk;
    }

    num_pat_lens = smpl_pat_lens[dir].size();
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

    std::vector<sidx_t> s_p(p + 1);
    s_p[0] = 0;
    s_p[p] = s;

    {
        std::array<sidx_t, (1 << 16) + 1> C2;
        C2[0] = 0;
        
        for (uint32_t c = 0; c < (1 << 16); c++) {
            if (SXIV2[dir][c].b == no_occ) {
                C2[c + 1] = C2[c];
            } else {
                C2[c + 1] = SXIV2[dir][c].e + 1;
            }
        }

        for (uint16_t i_p = 1; i_p < p; i_p++) {
            uint16_t c = bin_search_min_geq<sidx_t, uint16_t>(
                i_p * (s / p), 0, (1 << 16) - 1, [&](uint16_t c){return C2[c];});
            s_p[i_p] = C2[c];
        }
    }
    
    struct __attribute__((packed)) sxa_iv_fp_t {
        sxa_interval_t iv;
        uint64_t fp;
    };

    std::vector<std::vector<std::vector<sxa_iv_fp_t>>> SXIVX_vec(num_pat_lens,
        std::vector<std::vector<sxa_iv_fp_t>>(p));

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        sidx_t b = s_p[i_p];
        sidx_t e = s_p[i_p + 1];

        std::vector<sidx_t> sxivx_sizes(num_pat_lens, 0);

        for (sidx_t i = b + 1; i <= e; i++) {
            pos_t lcx = SLCX[i];

            for (int16_t j = 2; j < num_pat_lens; j++) {
                if (lcx < smpl_pat_lens[dir][j]) {
                    sxivx_sizes[j]++;
                }
            }
        }

        for (int16_t j = 2; j < num_pat_lens; j++) {
            SXIVX_vec[j][i_p].reserve(sxivx_sizes[j]);
        }
    }
    
    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        sidx_t b = s_p[i_p];
        sidx_t e = s_p[i_p + 1];
        std::vector<sidx_t> sxiv_b(num_pat_lens, b);

        for (sidx_t i = b + 1; i <= e; i++) {
            pos_t lcx = SLCX[i];

            for (int16_t j = num_pat_lens - 1; j >= 2; j--) {
                pos_t len = smpl_pat_lens[dir][j];
                if (lcx >= len) break;
                pos_t pos_im1 = S[SXA<dir>(i - 1)];

                if (is_pos_in_T<dir>(pos_im1, len - 1)) {
                    SXIVX_vec[j][i_p].emplace_back(sxa_iv_fp_t {
                        .iv = {.b = sxiv_b[j], .e = i - 1},
                        .fp = rks.substring<dir>(pos_im1, len)
                    });
                }

                sxiv_b[j] = i;
            }
        }
    }

    #pragma omp parallel for num_threads(p)
    for (int16_t j = 2; j < num_pat_lens; j++) {
        SXIVX<dir>()[j].reserve(pat_len_rank[j]);

        for (std::vector<sxa_iv_fp_t>& vec : SXIVX_vec[j]) {
            for (sxa_iv_fp_t& val : vec) {
                SXIVX<dir>()[j].emplace_with_hash(val.fp,
                    sxa_interval_t {.b = val.iv.b, .e = val.iv.e});
            }
            
            vec.clear();
            vec.shrink_to_fit();
        }
    }

    #ifndef NDEBUG
    std::vector<sidx_t> sxiv_b(num_pat_lens, 0);

    for (sidx_t i = 1; i <= s; i++) {
        pos_t lcx = SLCX[i];

        for (int16_t j = num_pat_lens - 1; j >= 2; j--) {
            pos_t len = smpl_pat_lens[dir][j];
            if (lcx >= len) break;
            pos_t pos_im1 = S[SXA<dir>(i - 1)];

            if (is_pos_in_T<dir>(pos_im1, len - 1)) {
                auto it = SXIVX<dir>()[j].find(pos_to_interval(pos_im1));
                assert(it != SXIVX<dir>()[j].end());
                assert(it->b == sxiv_b[j]);
                assert(it->e == i - 1);
            }
            
            sxiv_b[j] = i;
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