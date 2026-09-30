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

#include <lz77_sss/data_structures/sample_index/sample_index.hpp>

template <typename text_t, typename lce_r_t, typename array_t>
template <direction dir>
inline bool sample_index<text_t, lce_r_t, array_t>::cmp_sample_lex(uint64_t i, uint64_t j) const
{
    if (i == j) [[unlikely]] return false;
    return cmp_lex<dir>(S[i], S[j], lce<dir>(S[i], S[j], max_patt_len_left));
}

template <typename text_t, typename lce_r_t, typename array_t>
void sample_index<text_t, lce_r_t, array_t>::build_xiv_s_1_2(uint16_t p, bool log)
{
    auto time = now();

    if (log) {
        std::cout << "precomputing PA_C- and SA_C intervals for all 1- and 2-mers" << std::flush;
    }

    #pragma omp parallel for num_threads(p)
    for (uint16_t c = 0; c < 256; c++) {
        XIV_S_1[c].b = bin_search_min_geq<bool, uint64_t>(
            true, 0, s, [&](uint64_t i) {
                return i == s || T[S[SA_S[i]]] >= c;
            });

        if (XIV_S_1[c].b == s || T[S[SA_S[XIV_S_1[c].b]]] != c) {
            XIV_S_1[c].b = no_occ;
            XIV_S_1[c].e = no_occ;
        } else {
            XIV_S_1[c].e = bin_search_max_lt<bool, uint64_t>(
                true, XIV_S_1[c].b, s - 1, [&](uint64_t i) {
                    return T[S[SA_S[i]]] > c;
                });
        }
    }

    for (uint16_t c1 = 0; c1 < 256; c1++) {
        #pragma omp parallel for num_threads(p)
        for (uint16_t c2 = 0; c2 < 256; c2++) {
            uint16_t c_l = (c1 << 8) | c2;
            uint16_t c_r = (c2 << 8) | c1;

            if (XIV_S_1[c1].b == no_occ) {
                XIV_S_2[LEFT][c_l].b = no_occ;
                XIV_S_2[LEFT][c_l].e = no_occ;
                XIV_S_2[RIGHT][c_r].b = no_occ;
                XIV_S_2[RIGHT][c_r].e = no_occ;
            } else {
                XIV_S_2[LEFT][c_l].b = bin_search_min_geq<bool, uint64_t>(
                    true, XIV_S_1[c1].b, XIV_S_1[c1].e + 1, [&](uint64_t i) {
                        return i == XIV_S_1[c1].e + 1 || (S[PA_S[i]] != 0 && T[S[PA_S[i]] - 1] >= c2);
                    });

                if (XIV_S_2[LEFT][c_l].b == XIV_S_1[c1].e + 1 || T[S[PA_S[XIV_S_2[LEFT][c_l].b]] - 1] != c2) {
                    XIV_S_2[LEFT][c_l].b = no_occ;
                    XIV_S_2[LEFT][c_l].e = no_occ;
                } else {
                    XIV_S_2[LEFT][c_l].e = bin_search_max_lt<bool, uint64_t>(
                        true, XIV_S_2[LEFT][c_l].b, XIV_S_1[c1].e, [&](uint64_t i) {
                            return S[PA_S[i]] != 0 && T[S[PA_S[i]] - 1] > c2;
                        });
                }

                XIV_S_2[RIGHT][c_r].b = bin_search_min_geq<bool, uint64_t>(
                    true, XIV_S_1[c1].b, XIV_S_1[c1].e + 1, [&](uint64_t i) {
                        return i == XIV_S_1[c1].e + 1 || (S[SA_S[i]] != n - 1 && T[S[SA_S[i]] + 1] >= c2);
                    });

                if (XIV_S_2[RIGHT][c_r].b == XIV_S_1[c1].e + 1 || T[S[SA_S[XIV_S_2[RIGHT][c_r].b]] + 1] != c2) {
                    XIV_S_2[RIGHT][c_r].b = no_occ;
                    XIV_S_2[RIGHT][c_r].e = no_occ;
                } else {
                    XIV_S_2[RIGHT][c_r].e = bin_search_max_lt<bool, uint64_t>(
                        true, XIV_S_2[RIGHT][c_r].b, XIV_S_1[c1].e, [&](uint64_t i) {
                            return S[SA_S[i]] != n - 1 && T[S[SA_S[i]] + 1] > c2;
                        });
                }
            }
        }
    }

    if (log) {
        time = log_runtime(time);
    }
}

template <typename text_t, typename lce_r_t, typename array_t>
template <direction dir>
void sample_index<text_t, lce_r_t, array_t>::build_interval_hash_sets(uint64_t max_smpl_len, uint16_t p, bool log)
{
    auto time = now();
    uint64_t baseline = malloc_count_current();

    if (log) {
        std::cout << "building LC" << (dir == LEFT ? "S" : "P") << "_C" << std::flush;
    }

    const uint64_t lcx_cap = std::min<uint64_t>({ max_smpl_len, n, uint64_t(1) << 16 });
    bit_aligned_vector LCX_S(s + 1, lcx_cap);

    parallel_chunks(s + 1, p, [&](uint64_t beg, uint64_t end) {
        LCX_S.fill(beg, end, [&](uint64_t i) {
            if (i == 0 || i >= s) [[unlikely]] return uint64_t(0);

            return std::min<uint64_t>(lcx_cap, lce<dir>(
                S[XA_S<dir>(i - 1)], S[XA_S<dir>(i)], max_patt_len_left));
        });
    });

    std::vector<uint64_t> lcx_ranks(lcx_cap + 2, 0);
    uint64_t lcx_max = 0;

    {
        std::vector<std::vector<uint64_t>> lcx_cnt(p, std::vector<uint64_t>(lcx_cap + 1, 0));

        #pragma omp parallel for num_threads(p) schedule(dynamic, 65536)
        for (uint64_t i = 0; i < s; i++) {
            lcx_cnt[omp_get_thread_num()][LCX_S[i]]++;
        }

        uint64_t sum = 0;

        for (uint64_t v = 0; v <= lcx_cap; v++) {
            uint64_t cnt = 0;
            for (uint16_t i_p = 0; i_p < p; i_p++) cnt += lcx_cnt[i_p][v];
            if (cnt != 0) lcx_max = v;
            lcx_ranks[v] = sum;
            sum += cnt;
        }

        lcx_ranks[lcx_cap + 1] = sum;
    }

    auto lcx_rank_of = [&](uint64_t value) {
        return std::min<uint64_t>(lcx_ranks[std::min<uint64_t>(value, lcx_cap + 1)], s - 1);
    };

    auto lcx_value_at = [&](uint64_t rank) {
        return bin_search_min_geq<uint64_t, uint64_t>(rank + 1, 0, lcx_cap,
            [&](uint64_t v) { return lcx_ranks[v + 1]; });
    };

    if (log) {
        time = log_runtime(time);
    }

    max_smpl_len = std::min<uint64_t>(lcx_max, max_smpl_len);
    uint64_t lcx_s_rng_min = lcx_rank_of(3);
    uint64_t lcx_s_rng_max = lcx_rank_of(max_smpl_len);

    double max_num_ivs = 2.0 * s;
    double lcx_s_rng = lcx_s_rng_max - lcx_s_rng_min;
    smpl_patt_lens[dir] = { 1, 2 };
    XIV_S[dir].resize(2);
    if (byte_text && lcx_s_rng_min >= lcx_s_rng_max) return;
    uint64_t num_patt_lens = lcx_s_rng_min >= lcx_s_rng_max ? 2 : std::min<uint64_t>(max_smpl_len - 2,
        2 + std::floor((2.0 * max_num_ivs + lcx_s_rng_max - lcx_s_rng_min) /
        (double)(lcx_s_rng_max + lcx_s_rng_min)));
    std::vector<uint64_t> patt_len_ranks(num_patt_lens, 0);

    for (uint64_t i = 2; i < num_patt_lens; i++) {
        double rel_lcx_rank = (i - 1) / (double)(num_patt_lens - 2);
        uint64_t lcx_rank = std::floor((double)(lcx_s_rng_min) + rel_lcx_rank * lcx_s_rng);
        uint64_t len = std::max(lcx_value_at(lcx_rank), smpl_patt_lens[dir].back() + 1);
        if (len > max_smpl_len) break;

        smpl_patt_lens[dir].emplace_back(len);
        patt_len_ranks[i] = lcx_rank_of(len);
        XIV_S[dir].emplace_back(interval_hash_set());
    }

    uint64_t max_ivs_to_add = max_num_ivs * 0.2;

    if (num_patt_lens > 2 && patt_len_ranks[2] < max_ivs_to_add) {
        uint64_t added_ivs = 0;

        for (uint64_t len = 3; true; len++) {
            if (contains(smpl_patt_lens[dir], len)) continue;

            uint64_t lcx_rank = lcx_rank_of(len);
            if (len > max_smpl_len || added_ivs + lcx_rank > max_ivs_to_add) break;

            added_ivs += lcx_rank;
            smpl_patt_lens[dir].insert(smpl_patt_lens[dir].begin() + len - 1, len);
            patt_len_ranks.insert(patt_len_ranks.begin() + len - 1, lcx_rank);
            XIV_S[dir].insert(XIV_S[dir].begin() + len - 1, interval_hash_set());
        }
    }

    num_patt_lens = smpl_patt_lens[dir].size();
    lcx_ranks.clear();
    lcx_ranks.shrink_to_fit();
    if (num_patt_lens <= first_hashed_len_idx) return;

    if (log) {
        std::cout << "chose " << num_patt_lens - first_hashed_len_idx << " pattern"
            << " lengths in the range [" << smpl_patt_lens[dir][first_hashed_len_idx]
            << ", " << smpl_patt_lens[dir].back() << "]";
        time = log_runtime(time);
        std::cout << "sampling " << (dir == LEFT ? "P" : "S")
            << "A_C intervals" << std::flush;
    }

    const uint64_t len_min = smpl_patt_lens[dir][first_hashed_len_idx];
    const uint64_t num_chunks = std::max<uint64_t>(1, std::min<uint64_t>(s, uint64_t(p) * 16));
    const uint64_t chunk_target = div_ceil(s, num_chunks);
    std::vector<uint64_t> sect_beg;
    sect_beg.reserve(num_chunks + 1);
    sect_beg.emplace_back(0);

    for (uint64_t i = 1; i < s; i++) {
        if (i - sect_beg.back() >= chunk_target && LCX_S[i] < len_min) sect_beg.emplace_back(i);
    }

    sect_beg.emplace_back(s);
    const uint64_t num_sects = sect_beg.size() - 1;

    std::vector<std::vector<uint64_t>> xiv_s_cnt(p, std::vector<uint64_t>(num_patt_lens, 0));
    std::vector<std::vector<uint64_t>> xiv_s_max(p, std::vector<uint64_t>(num_patt_lens, 0));
    const uint64_t len_max = smpl_patt_lens[dir].back();

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t sect = 0; sect < num_sects; sect++) {
        const uint16_t i_p = omp_get_thread_num();
        const uint64_t beg = sect_beg[sect];
        const uint64_t end = sect_beg[sect + 1];
        std::vector<uint64_t> iv_beg(num_patt_lens, beg);
        std::vector<uint64_t>& cnt = xiv_s_cnt[i_p];
        std::vector<uint64_t>& max_width = xiv_s_max[i_p];

        for (uint64_t i = beg + 1; i <= end; i++) {
            const uint64_t lcx = LCX_S[i];
            if (lcx >= len_max) continue;
            const uint64_t pos_im1 = S[XA_S<dir>(i - 1)];

            for (int64_t j = num_patt_lens - 1; j >= int64_t(first_hashed_len_idx); j--) {
                const uint64_t len = smpl_patt_lens[dir][j];
                if (lcx >= len) break;

                if (is_pos_in_T<dir>(pos_im1, len - 1)) {
                    cnt[j]++;
                    max_width[j] = std::max<uint64_t>(max_width[j], i - 1 - iv_beg[j]);
                }

                iv_beg[j] = i;
            }
        }
    }

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t j = first_hashed_len_idx; j < num_patt_lens; j++) {
        uint64_t cnt = 0;
        uint64_t max_width = 0;

        for (uint16_t i_p = 0; i_p < p; i_p++) {
            cnt += xiv_s_cnt[i_p][j];
            max_width = std::max<uint64_t>(max_width, xiv_s_max[i_p][j]);
        }

        XIV_S[dir][j] = interval_hash_set(cnt, s, max_width);
    }

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t sect = 0; sect < num_sects; sect++) {
        const uint64_t beg = sect_beg[sect];
        const uint64_t end = sect_beg[sect + 1];
        std::vector<uint64_t> iv_beg(num_patt_lens, beg);

        for (uint64_t i = beg + 1; i <= end; i++) {
            const uint64_t lcx = LCX_S[i];
            if (lcx >= len_max) continue;
            const uint64_t pos_im1 = S[XA_S<dir>(i - 1)];

            for (int64_t j = num_patt_lens - 1; j >= int64_t(first_hashed_len_idx); j--) {
                const uint64_t len = smpl_patt_lens[dir][j];
                if (lcx >= len) break;

                if (is_pos_in_T<dir>(pos_im1, len - 1)) {
                    XIV_S[dir][j].insert_parallel(
                        RKS.template substr_fp<dir>(pos_im1, len), iv_beg[j], i - 1);
                }

                iv_beg[j] = i;
            }
        }
    }

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t j = first_hashed_len_idx; j < num_patt_lens; j++) {
        XIV_S[dir][j].finish();
    }

    if (log) {
        std::string phase = dir == LEFT ? "pa_c_samples" : "sa_c_samples";
        record_phase_time(phase, time_diff_ns(time, now()));
        std::cout << " (" << format_size(malloc_count_current() - baseline) << ")";
        time = log_runtime(time);
    }
}

template <typename text_t, typename lce_r_t, typename array_t>
void sample_index<text_t, lce_r_t, array_t>::build(
    const text_t& T,
    uint64_t n,
    const array_t& S,
    const lce_r_t& LCE_R,
    interval_samples mode,
    uint64_t rks_sample_rate,
    uint16_t p,
    bool log,
    uint64_t max_patt_len_left,
    uint64_t max_smpl_len_right)
{
    uint64_t baseline_bytes = malloc_count_current();
    auto time = now();

    this->T = T;
    this->S = array_view(S);
    this->LCE_R = &LCE_R;
    this->max_patt_len_left = max_patt_len_left;
    this->n = n;
    s = S.size();

    #ifndef NDEBUG
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 1; i < s; i++) {
        assert(S[i] > S[i - 1]);
    }
    #endif

    if (log) {
        std::cout << "building PA_C" << std::flush;
    }

    {
        std::vector<uint40_t> sorted;
        no_init_resize(sorted, s);

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < s; i++) {
            sorted[i] = i;
        }

        ips4o::parallel::sort(sorted.begin(), sorted.end(),
            [&](uint64_t i, uint64_t j) {
                return cmp_sample_lex<LEFT>(i, j);
            });

        PA_S = bit_aligned_vector(s, s);

        parallel_chunks(s, p, [&](uint64_t beg, uint64_t end) {
            PA_S.fill(beg, end, [&](uint64_t i) { return uint64_t(sorted[i]); });
        });
    }

    #ifndef NDEBUG
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 1; i < s; i++) {
        assert(!cmp_sample_lex<LEFT>(PA_S[i], PA_S[i - 1]));
    }
    #endif

    if (log) {
        record_phase_time("pa_c", time_diff_ns(time, now()));
        std::cout << " (" << format_size(PA_S.size_in_bytes()) << ")" << std::flush;
        time = log_runtime(time);
        std::cout << "building SA_C" << std::flush;
    }

    {
        std::vector<uint40_t> sorted;
        no_init_resize(sorted, s);

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < s; i++) {
            sorted[i] = i;
        }

        if constexpr (requires { this->LCE_R->get_isa_s(); this->LCE_R->get_sync_set(); this->LCE_R->tau(); }) {
            const auto& sync = LCE_R.get_sync_set();
            const auto& isa = LCE_R.get_isa_s();
            const uint64_t num_sync = sync.size();
            const uint64_t max_off = 3 * LCE_R.tau();
            constexpr uint32_t no_off = std::numeric_limits<uint32_t>::max();
            std::vector<uint32_t> off;
            std::vector<uint40_t> rank;
            no_init_resize(off, s);
            no_init_resize(rank, s);

            parallel_chunks(s, p, [&](uint64_t beg, uint64_t end) {
                uint64_t j = bin_search_min_geq<uint64_t, uint64_t>(
                    S[beg], 0, num_sync, [&](uint64_t x) { return x == num_sync ? n : uint64_t(sync[x]); });

                for (uint64_t i = beg; i < end; i++) {
                    const uint64_t pos = S[i];
                    while (j < num_sync && sync[j] < pos) j++;

                    if (j == num_sync || sync[j] - pos > max_off) {
                        off[i] = no_off;
                    } else {
                        off[i] = sync[j] - pos;
                        rank[i] = isa[j];
                    }
                }
            });

            ips4o::parallel::sort(sorted.begin(), sorted.end(), [&](uint64_t i, uint64_t j) {
                if (i == j) [[unlikely]] return false;
                const uint32_t oi = off[i];

                if (oi != no_off && oi == off[j]) {
                    const uint64_t a = S[i];
                    const uint64_t b = S[j];
                    const uint64_t l = T.lce(a, b, oi);
                    if (l < oi) return T[a + l] < T[b + l];
                    return uint64_t(rank[i]) < uint64_t(rank[j]);
                }

                return cmp_sample_lex<RIGHT>(i, j);
            });
        } else {
            ips4o::parallel::sort(sorted.begin(), sorted.end(),
                [&](uint64_t i, uint64_t j) {
                    return cmp_sample_lex<RIGHT>(i, j);
                });
        }

        SA_S = bit_aligned_vector(s, s);

        parallel_chunks(s, p, [&](uint64_t beg, uint64_t end) {
            SA_S.fill(beg, end, [&](uint64_t i) { return uint64_t(sorted[i]); });
        });
    }

    #ifndef NDEBUG
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 1; i < s; i++) {
        assert(!cmp_sample_lex<RIGHT>(SA_S[i], SA_S[i - 1]));
    }
    #endif

    if (log) {
        record_phase_time("sa_c", time_diff_ns(time, now()));
        std::cout << " (" << format_size(SA_S.size_in_bytes()) << ")" << std::flush;
        time = log_runtime(time);
    }

    if (mode == interval_samples::use) {
        if (log) {
            std::cout << "building RKS" << std::flush;
        }

        RKS = rks_t(T, n, rks_sample_rate, 0, p);

        if (log) {
            record_phase_time("rks", time_diff_ns(time, now()));
            std::cout << " (" << format_size(RKS.size_in_bytes()) << ")" << std::flush;
            time = log_runtime(time);
        }

        if constexpr (byte_text) build_xiv_s_1_2(p, log);
        build_interval_hash_sets<LEFT>(max_patt_len_left, p, log);
        build_interval_hash_sets<RIGHT>(max_smpl_len_right, p, log);
    }

    byte_size = malloc_count_current() - baseline_bytes;
}
