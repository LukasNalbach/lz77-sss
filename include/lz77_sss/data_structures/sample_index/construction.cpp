#pragma once

#include <lz77_sss/data_structures/sample_index/sample_index.hpp>

template <typename pos_t, typename sidx_t, typename char_t, typename lce_r_t>
void sample_index<pos_t, sidx_t, char_t, lce_r_t>::build_xa_s_1_2_intervals(uint16_t p, bool log)
{
    auto time = now();

    if (log) {
        std::cout << "precomputing SPA1/2 and SSA1/2 intervals" << std::flush;
    }

    #pragma omp parallel for num_threads(p)
    for (uint16_t c = 0; c < 256; c++) {
        SIV_S_1[c].b = bin_search_min_geq<bool, sidx_t>(
            true, 0, s, [&](sidx_t i) {
                return i == s || char_to_uchar(T[S[SA_S[i]]]) >= c;
            });

        if (SIV_S_1[c].b == s || char_to_uchar(T[S[SA_S[SIV_S_1[c].b]]]) != c) {
            SIV_S_1[c].b = no_occ;
            SIV_S_1[c].e = no_occ;
        } else {
            SIV_S_1[c].e = bin_search_max_lt<bool, sidx_t>(
                true, SIV_S_1[c].b, s - 1, [&](sidx_t i) {
                    return char_to_uchar(T[S[SA_S[i]]]) > c;
                });
        }
    }

    for (uint16_t c1 = 0; c1 < 256; c1++) {
        #pragma omp parallel for num_threads(p)
        for (uint16_t c2 = 0; c2 < 256; c2++) {
            uint16_t c_l = (c1 << 8) | c2;
            uint16_t c_r = (c2 << 8) | c1;

            if (SIV_S_1[c1].b == no_occ) {
                XIV_S_2[LEFT][c_l].b = no_occ;
                XIV_S_2[LEFT][c_l].e = no_occ;
                XIV_S_2[RIGHT][c_r].b = no_occ;
                XIV_S_2[RIGHT][c_r].e = no_occ;
            } else {
                XIV_S_2[LEFT][c_l].b = bin_search_min_geq<bool, sidx_t>(
                    true, SIV_S_1[c1].b, SIV_S_1[c1].e + 1, [&](sidx_t i) {
                        return i == SIV_S_1[c1].e + 1 || (S[PA_S[i]] != 0 && char_to_uchar(T[S[PA_S[i]] - 1]) >= c2);
                    });

                if (XIV_S_2[LEFT][c_l].b == SIV_S_1[c1].e + 1 || char_to_uchar(T[S[PA_S[XIV_S_2[LEFT][c_l].b]] - 1]) != c2) {
                    XIV_S_2[LEFT][c_l].b = no_occ;
                    XIV_S_2[LEFT][c_l].e = no_occ;
                } else {
                    XIV_S_2[LEFT][c_l].e = bin_search_max_lt<bool, sidx_t>(
                        true, XIV_S_2[LEFT][c_l].b, SIV_S_1[c1].e, [&](sidx_t i) {
                            return S[PA_S[i]] != 0 && char_to_uchar(T[S[PA_S[i]] - 1]) > c2;
                        });
                }

                XIV_S_2[RIGHT][c_r].b = bin_search_min_geq<bool, sidx_t>(
                    true, SIV_S_1[c1].b, SIV_S_1[c1].e + 1, [&](sidx_t i) {
                        return i == SIV_S_1[c1].e + 1 || (S[SA_S[i]] != n - 1 && char_to_uchar(T[S[SA_S[i]] + 1]) >= c2);
                    });

                if (XIV_S_2[RIGHT][c_r].b == SIV_S_1[c1].e + 1 || char_to_uchar(T[S[SA_S[XIV_S_2[RIGHT][c_r].b]] + 1]) != c2) {
                    XIV_S_2[RIGHT][c_r].b = no_occ;
                    XIV_S_2[RIGHT][c_r].e = no_occ;
                } else {
                    XIV_S_2[RIGHT][c_r].e = bin_search_max_lt<bool, sidx_t>(
                        true, XIV_S_2[RIGHT][c_r].b, SIV_S_1[c1].e, [&](sidx_t i) {
                            return S[SA_S[i]] != n - 1 && char_to_uchar(T[S[SA_S[i]] + 1]) > c2;
                        });
                }
            }
        }
    }

    if (log) {
        time = log_runtime(time);
    }
}

template <typename pos_t, typename sidx_t, typename char_t, typename lce_r_t>
template <direction dir>
void sample_index<pos_t, sidx_t, char_t, lce_r_t>::build_samples(pos_t max_smpl_len, uint16_t p, bool log)
{
    auto time = now();

    if (log) {
        std::cout << "building LC" << (dir == LEFT ? "S" : "P") << "_S" << std::flush;
    }

    std::vector<pos_t> LCX_S;
    no_init_resize(LCX_S, s + 1);
    LCX_S[0] = 0;
    LCX_S[s] = 0;

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 1; i < s; i++) {
        LCX_S[i] = lce<dir>(
            S[XA_S<dir>(i - 1)],
            S[XA_S<dir>(i)],
            max_patt_len_left);
    }

    if (log) {
        time = log_runtime(time);
        std::cout << "sorting LC" << (dir == LEFT ? "S" : "P") << "_S" << std::flush;
    }

    uint64_t baseline = malloc_count_current();
    std::vector<pos_t> LCX_S_sorted = LCX_S;
    ips4o::parallel::sort(LCX_S_sorted.begin(), LCX_S_sorted.end() - 1);

    if (log) {
        time = log_runtime(time);
    }

    max_smpl_len = std::min<pos_t>(LCX_S_sorted[s - 1], max_smpl_len);

    sidx_t lcx_s_rng_min = bin_search_min_geq<pos_t, sidx_t>(3,
        0, s - 1, [&](sidx_t i) { return LCX_S_sorted[i]; });
    sidx_t lcx_s_rng_max = bin_search_min_geq<pos_t, sidx_t>(max_smpl_len,
        0, s - 1, [&](sidx_t i) { return LCX_S_sorted[i]; });

    double max_num_samples = 2.0 * s;
    double lcx_s_rng = lcx_s_rng_max - lcx_s_rng_min;
    smpl_pat_lens[dir] = { 1, 2 };
    XIV_S<dir>().resize(2);
    if (lcx_s_rng_min >= lcx_s_rng_max) return;
    uint64_t num_pat_lens = std::min<uint64_t>(max_smpl_len - 2,
        2 + std::floor((2.0 * max_num_samples) / double { lcx_s_rng_min + lcx_s_rng_max }));
    std::vector<sidx_t> pat_len_rank(num_pat_lens, 0);

    for (uint64_t i = 2; i < num_pat_lens; i++) {
        double rel_lcx_s_rnk = (i - 1) / double { num_pat_lens - 2 };
        pos_t lcx_s_rnk = std::floor(double { lcx_s_rng_min } + rel_lcx_s_rnk * lcx_s_rng);
        pos_t len = std::max(LCX_S_sorted[lcx_s_rnk], smpl_pat_lens[dir].back() + 1);
        if (len > max_smpl_len) break;

        smpl_pat_lens[dir].emplace_back(len);
        lcx_s_rnk = bin_search_min_geq<pos_t, sidx_t>(len,
            0, s - 1, [&](sidx_t i) {return LCX_S_sorted[i];});

        XIV_S<dir>().emplace_back(xiv_s_map_t<dir>(0,
            xiv_s_hash<dir>(this, len),
            xiv_s_eq<dir>(this, len)));

        pat_len_rank[i] = lcx_s_rnk;
    }

    num_pat_lens = smpl_pat_lens[dir].size();
    LCX_S_sorted.clear();
    LCX_S_sorted.shrink_to_fit();

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
            if (XIV_S_2[dir][c].b == no_occ) {
                C2[c + 1] = C2[c];
            } else {
                C2[c + 1] = XIV_S_2[dir][c].e + 1;
            }
        }

        for (uint16_t i_p = 1; i_p < p; i_p++) {
            uint16_t c = bin_search_min_geq<sidx_t, uint16_t>(i_p * (s / p),
                0, (1 << 16) - 1, [&](uint16_t c) { return C2[c]; });
            s_p[i_p] = C2[c];
        }
    }

    struct __attribute__((packed)) sxa_iv_fp_t {
        interval_t iv;
        uint32_t fp;
    };

    std::vector<std::vector<std::vector<sxa_iv_fp_t>>> XIV_S_vec(num_pat_lens,
        std::vector<std::vector<sxa_iv_fp_t>>(p));

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        sidx_t b = s_p[i_p];
        sidx_t e = s_p[i_p + 1];

        std::vector<sidx_t> xiv_s_sizes(num_pat_lens, 0);

        for (sidx_t i = b + 1; i <= e; i++) {
            pos_t lcx = LCX_S[i];

            for (uint64_t j = 2; j < num_pat_lens; j++) {
                if (lcx < smpl_pat_lens[dir][j]) {
                    xiv_s_sizes[j]++;
                }
            }
        }

        for (uint64_t j = 2; j < num_pat_lens; j++) {
            XIV_S_vec[j][i_p].reserve(xiv_s_sizes[j]);
        }
    }

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        sidx_t b = s_p[i_p];
        sidx_t e = s_p[i_p + 1];
        std::vector<sidx_t> xiv_s_b(num_pat_lens, b);

        for (sidx_t i = b + 1; i <= e; i++) {
            pos_t lcx = LCX_S[i];

            for (int16_t j = num_pat_lens - 1; j >= 2; j--) {
                pos_t len = smpl_pat_lens[dir][j];
                if (lcx >= len) break;
                pos_t pos_im1 = S[XA_S<dir>(i - 1)];

                if (is_pos_in_T<dir>(pos_im1, len - 1)) {
                    XIV_S_vec[j][i_p].emplace_back(sxa_iv_fp_t {
                        .iv = { .b = xiv_s_b[j], .e = i - 1 },
                        .fp = rks.template substring<dir>(pos_im1, len) });
                }

                xiv_s_b[j] = i;
            }
        }
    }

    #pragma omp parallel for num_threads(p)
    for (uint64_t j = 2; j < num_pat_lens; j++) {
        XIV_S<dir>()[j].reserve(pat_len_rank[j]);

        for (std::vector<sxa_iv_fp_t>& vec : XIV_S_vec[j]) {
            for (sxa_iv_fp_t& val : vec) {
                XIV_S<dir>()[j].emplace_with_hash(val.fp,
                    interval_t { .b = val.iv.b, .e = val.iv.e });
            }

            vec.clear();
            vec.shrink_to_fit();
        }
    }

    #ifndef NDEBUG
    std::vector<sidx_t> xiv_s_b(num_pat_lens, 0);

    for (sidx_t i = 1; i <= s; i++) {
        pos_t lcx = LCX_S[i];

        for (uint64_t j = num_pat_lens - 1; j >= 2; j--) {
            pos_t len = smpl_pat_lens[dir][j];
            if (lcx >= len) break;
            pos_t pos_im1 = S[XA_S<dir>(i - 1)];

            if (is_pos_in_T<dir>(pos_im1, len - 1)) {
                auto it = XIV_S<dir>()[j].find(pos_to_interval(pos_im1));
                assert(it != XIV_S<dir>()[j].end());
                assert(it->b == xiv_s_b[j]);
                assert(it->e == i - 1);
            }

            xiv_s_b[j] = i;
        }

        if (lcx < 2) {
            pos_t pos_im1 = S[XA_S<dir>(i - 1)];

            if (is_pos_in_T<dir>(pos_im1, 1)) {
                uint16_t pat = val_offs<uint16_t, dir, 0>(T, pos_im1);
                assert(XIV_S_2[dir][pat].b == xiv_s_b[1]);
                assert(XIV_S_2[dir][pat].e == i - 1);
            }

            xiv_s_b[1] = i;

            if (dir == LEFT && lcx < 1) {
                uint8_t pat = val_offs<uint8_t, RIGHT, 0>(T, pos_im1);
                assert(SIV_S_1[pat].b == xiv_s_b[0]);
                assert(SIV_S_1[pat].e == i - 1);
                xiv_s_b[0] = i;
            }
        }
    }
    #endif

    if (log) {
        std::string phase = dir == LEFT ? "spa_samples" : "ssa_samples";
        log_phase(phase, time_diff_ns(time, now()));
        std::cout << " (" << format_size(malloc_count_current() - baseline) << ")";
        time = log_runtime(time);
    }
}