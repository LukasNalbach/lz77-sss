#pragma once

template <typename pos_t, typename sidx_t, typename lce_r_t>
void sample_index<pos_t, sidx_t, lce_r_t>::build_spa12_intervals(bool log) {
    auto time = now();

    if (log) {
        std::cout << "precomputing SPA1/2 and SSA1/2 intervals" << std::flush;
    }

    for (uint16_t c = 0; c < 256; c++) {
        SCIV[c].b = bin_search_min_geq<bool, sidx_t>(
            true, 0, s, [&](sidx_t i){
                return i == s || T[S[SSA[i]]] >= c;
        });

        if (SCIV[c].b == s || T[S[SSA[SCIV[c].b]]] != c) {
            SCIV[c].b = no_occ;
            SCIV[c].e = no_occ;
        } else {
            SCIV[c].e = bin_search_max_lt<bool, sidx_t>(
                true, SCIV[c].b, s - 1, [&](sidx_t i){
                    return T[S[SSA[i]]] > c;
            });

            #ifndef NDEBUG
            sidx_t b = SCIV[c].b;
            sidx_t e = SCIV[c].e;

            assert(b == 0     || T[S[SPA[b - 1]]] < c);
            assert(e == s - 1 || T[S[SPA[e + 1]]] > c);

            assert(b == 0     || T[S[SSA[b - 1]]] < c);
            assert(e == s - 1 || T[S[SSA[e + 1]]] > c);

            for (sidx_t i = b; i <= e; i++) {
                assert(T[S[SSA[i]]] == c);
                assert(T[S[SPA[i]]] == c);
            }
            #endif
        }
    }

    for (uint16_t c1 = 0; c1 < 256; c1++) {
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
                        return i == SCIV[c1].e + 1 || (S[SPA[i]] != 0 && T[S[SPA[i]] - 1] >= c2);
                });

                if (SXIV2[LEFT][c_l].b == SCIV[c1].e + 1 || T[S[SPA[SXIV2[LEFT][c_l].b]] - 1] != c2) {
                    SXIV2[LEFT][c_l].b = no_occ;
                    SXIV2[LEFT][c_l].e = no_occ;
                } else {
                    SXIV2[LEFT][c_l].e = bin_search_max_lt<bool, sidx_t>(
                        true, SXIV2[LEFT][c_l].b, SCIV[c1].e, [&](sidx_t i){
                            return S[SPA[i]] != 0 && T[S[SPA[i]] - 1] > c2;
                    });

                    #ifndef NDEBUG
                    sidx_t b_c1 = SCIV[c1].b;
                    sidx_t e_c1 = SCIV[c1].e;

                    sidx_t b = SXIV2[LEFT][c_l].b;
                    sidx_t e = SXIV2[LEFT][c_l].e;

                    assert(b >= b_c1);
                    assert(e <= e_c1);

                    assert(b == 0     || b == b_c1 || S[SPA[b - 1]] == 0 || T[S[SPA[b - 1]] - 1] < c2);
                    assert(e == s - 1 || e == e_c1 || S[SPA[e + 1]] == 0 || T[S[SPA[e + 1]] - 1] > c2);

                    for (sidx_t i = b; i <= e; i++) {
                        assert(T[S[SPA[i]] - 1] == c2);
                        assert(T[S[SPA[i]]] == c1);
                        assert(*reinterpret_cast<const uint16_t*>(&T[S[SPA[i]] - 1]) == c_l);
                    }
                    #endif
                }

                SXIV2[RIGHT][c_r].b = bin_search_min_geq<bool, sidx_t>(
                    true, SCIV[c1].b, SCIV[c1].e + 1, [&](sidx_t i){
                        return i == SCIV[c1].e + 1 || (S[SSA[i]] != n - 1 && T[S[SSA[i]] + 1] >= c2);
                });

                if (SXIV2[RIGHT][c_r].b == SCIV[c1].e + 1 || T[S[SSA[SXIV2[RIGHT][c_r].b]] + 1] != c2) {
                    SXIV2[RIGHT][c_r].b = no_occ;
                    SXIV2[RIGHT][c_r].e = no_occ;
                } else {
                    SXIV2[RIGHT][c_r].e = bin_search_max_lt<bool, sidx_t>(
                        true, SXIV2[RIGHT][c_r].b, SCIV[c1].e, [&](sidx_t i){
                            return S[SSA[i]] != n - 1 && T[S[SSA[i]] + 1] > c2;
                    });

                    #ifndef NDEBUG
                    sidx_t b_c1 = SCIV[c1].b;
                    sidx_t e_c1 = SCIV[c1].e;

                    sidx_t b = SXIV2[RIGHT][c_r].b;
                    sidx_t e = SXIV2[RIGHT][c_r].e;

                    assert(b >= b_c1);
                    assert(e <= e_c1);

                    assert(b == 0     || b == b_c1 || S[SSA[b - 1]] == n - 1 || T[S[SSA[b - 1]] + 1] < c2);
                    assert(e == s - 1 || e == e_c1 || S[SSA[e + 1]] == n - 1 || T[S[SSA[e + 1]] + 1] > c2);

                    for (sidx_t i = b; i <= e; i++) {
                        assert(T[S[SSA[i]]] == c1);
                        assert(T[S[SSA[i]] + 1] == c2);
                        assert(*reinterpret_cast<const uint16_t*>(&T[S[SSA[i]]]) == c_r);
                    }
                    #endif
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
void sample_index<pos_t, sidx_t, lce_r_t>::build_samples(pos_t typ_lce_r, bool log) {
    auto time = now();

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

    if (SLCX_sorted[s - 1] >= 3) {
        if (log) {
            time = log_runtime(time);
        }

        sidx_t slcx_rng_min = bin_search_min_geq<pos_t, sidx_t>(
            3, 0, s - 1, [&](sidx_t i){
                return SLCX_sorted[i];
        });

        pos_t max_sampled_len = std::min<pos_t>(
            SLCX_sorted[s - 1],
            dir == LEFT ? max_lce_l : typ_lce_r
        );

        sidx_t slcx_rng_max = bin_search_min_geq<pos_t, sidx_t>(
            max_sampled_len, 0, s - 1, [&](sidx_t i){
                return SLCX_sorted[i];
        });

        std::vector<pos_t> pat_lens;
        double max_num_samples = (dir == LEFT ? 4 : 4) * s;
        double slcx_rng = slcx_rng_max - slcx_rng_min;
        int16_t num_pat_lens = std::max<int16_t>(1, std::floor(
            (max_num_samples - slcx_rng / 2.0) /
            (double{slcx_rng_min} + slcx_rng / 2.0)
        ));

        for (int16_t i = 0; i < num_pat_lens; i++) {
            double rel_slcx_rnk = (i + 1.0) / double{num_pat_lens};
            pos_t slcx_rnk = slcx_rng_min + rel_slcx_rnk * slcx_rng;
            pos_t len = SLCX_sorted[slcx_rnk];

            if (!pat_lens.empty()) {
                while (len <= pat_lens.back()) {
                    len++;
                }

                if (len > max_sampled_len) {
                    break;
                }
            }

            pat_lens.emplace_back(len);
            slcx_rnk = bin_search_min_geq<pos_t, sidx_t>(
                len, 0, s - 1, [&](sidx_t i){
                    return SLCX_sorted[i];
            });

            SXIVX<dir>().emplace_back(sxivx_map_t<dir>(0,
                sxivx_hash<dir>(T, rks, len),
                sxivx_eq<dir>(*this, len)
            ));

            SXIVX<dir>().back().reserve(slcx_rnk);
        }

        num_pat_lens = pat_lens.size();
        smpl_pat_lens[dir].insert(smpl_pat_lens[dir].end(),
            pat_lens.begin(), pat_lens.end());
        SLCX_sorted.clear();
        SLCX_sorted.shrink_to_fit();
        std::vector<sidx_t> sxiv_b(num_pat_lens, 0);

        if (log) {
            std::cout << "chose " << num_pat_lens << " pattern"
                << " lengths in the range [" << pat_lens.front()
                << ", " << pat_lens.back() << "]";
            time = log_runtime(time);
            std::cout << "sampling S" << (dir == LEFT ? "P" : "S")
                << "A intervals" << std::flush;
        }

        for (sidx_t i = 1; i <= s; i++) {
            pos_t lcx = SLCX[i];

            for (int16_t idx = num_pat_lens - 1; idx >= 0; idx--) {
                pos_t len = pat_lens[idx];

                if (lcx >= len) {
                    break;
                }

                pos_t pos_im1 = S[SXA<dir>(i - 1)];

                if (is_pos_in_T<dir>(pos_im1, len - 1)) {
                    SXIVX<dir>()[idx].emplace(pos_im1,
                        sxa_interval_t{sxiv_b[idx], i - 1});
                }

                sxiv_b[idx] = i;
            }
        }

        #ifndef NDEBUG
        for (uint8_t idx = 0; idx < pat_lens.size(); idx++) {
            sxiv_b[idx] = 0;
        }

        for (sidx_t i = 1; i <= s; i++) {
            pos_t lcx = SLCX[i];

            for (int16_t idx = num_pat_lens - 1; idx >= 0; idx--) {
                pos_t len = pat_lens[idx];

                if (lcx >= len) {
                    break;
                }

                pos_t pos_im1 = S[SXA<dir>(i - 1)];

                if (is_pos_in_T<dir>(pos_im1, len - 1)) {
                    auto it = SXIVX<dir>()[idx].find(pos_im1);
                    assert(it != SXIVX<dir>()[idx].end());
                    assert(it.value().b == sxiv_b[idx]);
                    assert(it.value().e == i - 1);
                }
                
                sxiv_b[idx] = i;
            }
        }
        #endif

        if (log) {
            std::cout << " (" << format_size(malloc_count_current() - baseline) << ")";
            time = log_runtime(time);
        }
    }
}