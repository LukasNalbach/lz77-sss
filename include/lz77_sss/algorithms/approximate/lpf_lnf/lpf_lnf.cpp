#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <phrase_mode phr_mode>
void lz77_sss<pos_t>::factorizer<tau, char_t>::build_LPF_all()
{
    build_PSV_NSV_S();

    if (log) {
        std::cout << "building LPF" << std::flush;
    }

    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SA_S = LCE.get_sa_s();
    const std::vector<uint32_t>& ISA_S = LCE.get_isa_s();
    pos_t s = S.size();

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b = i_p * (n / p);
        pos_t e = i_p == p - 1 ? n : ((i_p + 1) * (n / p));

        uint32_t i_min = bin_search_min_geq<pos_t, uint32_t>(
            true, 0, s, [&](uint32_t i) { return i == s || S[i] >= b; });
        uint32_t i_max = bin_search_min_geq<pos_t, uint32_t>(
            true, 0, s, [&](uint32_t i) { return i == s || S[i] >= e; });

        lpf lst_sm_phr = {n, n, n};
        lpf lst_gr_phr = {n, n, n};

        for (uint32_t i = i_min; i < i_max; i++) {
            if (PSV_S[ISA_S[i]] != s) [[likely]] {
                pos_t beg = S[i];
                pos_t src = S[SA_S[PSV_S[ISA_S[i]]]];

                if (!(beg < lst_sm_phr.end &&
                    beg - src == lst_sm_phr.beg - lst_sm_phr.src
                )) {
                    pos_t end = S[i] + LCE_R(src, S[i]);

                    if (phr_mode == lpf_lnf_opt && src != 0 && S[i] != 0) {
                        pos_t lce_l = LCE_L(src - 1, S[i] - 1);
                        beg -= lce_l;
                        src -= lce_l;
                    }

                    #ifndef NDEBUG
                    assert(src < beg);

                    for (pos_t j = 0; j < end - beg; j++) {
                        assert(T[src + j] == T[beg + j]);
                    }
                    #endif

                    if (end - beg > 1) [[likely]] {
                        lst_sm_phr = { beg, end, src };
                        LPF[i_p].emplace_back(lst_sm_phr);
                    }
                }
            }

            if (NSV_S[ISA_S[i]] != s) [[likely]] {
                pos_t beg = S[i];
                pos_t src = S[SA_S[NSV_S[ISA_S[i]]]];

                if (!(beg < lst_gr_phr.end &&
                    beg - src == lst_gr_phr.beg - lst_gr_phr.src
                )) {
                    pos_t end = S[i] + LCE_R(src, S[i]);

                    if (phr_mode == lpf_lnf_opt && src != 0 && S[i] != 0) {
                        pos_t lce_l = LCE_L(src - 1, S[i] - 1);
                        beg -= lce_l;
                        src -= lce_l;
                    }

                    #ifndef NDEBUG
                    assert(src < beg);

                    for (pos_t j = 0; j < end - beg; j++) {
                        assert(T[src + j] == T[beg + j]);
                    }
                    #endif

                    if (end - beg > 1) [[likely]] {
                        lst_gr_phr = { beg, end, src };
                        LPF[i_p].emplace_back(lst_gr_phr);
                    }
                }
            }
        }
    }

    NSV_S.clear();
    PSV_S.clear();
    NSV_S.shrink_to_fit();
    PSV_S.shrink_to_fit();

    if (log) {
        log_phase("lpf", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <phrase_mode phr_mode>
void lz77_sss<pos_t>::factorizer<tau, char_t>::build_LNF_all()
{
    build_PGV_NGV_S();

    if (log) {
        std::cout << "building LNF" << std::flush;
    }

    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SA_S = LCE.get_sa_s();
    const std::vector<uint32_t>& ISA_S = LCE.get_isa_s();
    pos_t s = S.size();

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b = i_p == p - 1 ? 0 : (n - (i_p + 1) * (n / p));
        pos_t e = n - i_p * (n / p);

        uint32_t i_min = bin_search_min_geq<pos_t, uint32_t>(
            true, 0, s, [&](uint32_t i) { return i == s || S[i] >= b; });
        uint32_t i_max = bin_search_min_geq<pos_t, uint32_t>(
            true, 0, s, [&](uint32_t i) { return i == s || S[i] >= e; });
            
        lpf lst_sm_phr = {n, n, n};
        lpf lst_gr_phr = {n, n, n};

        for (uint32_t i = i_min; i < i_max; i++) {
            if (PGV_S[ISA_S[i]] != s) [[likely]] {
                pos_t src = S[SA_S[PGV_S[ISA_S[i]]]];
                pos_t beg = S[i];

                if (!(beg < lst_sm_phr.end &&
                    src - beg == lst_sm_phr.src - lst_sm_phr.beg
                )) {
                    pos_t end = S[i] + LCE_R(S[i], src);

                    if (phr_mode == lpf_lnf_opt && src != 0 && S[i] != 0) {
                        pos_t lce_l = LCE_L(src - 1, S[i] - 1);
                        beg -= lce_l;
                        src -= lce_l;
                    }

                    #ifndef NDEBUG
                    assert(src > beg);

                    for (pos_t j = 0; j < end - beg; j++) {
                        assert(T[src + j] == T[beg + j]);
                    }
                    #endif

                    if (end - beg > 1) [[likely]] {
                        lst_sm_phr = { beg, end, src };

                        LPF[i_p].emplace_back(lpf {
                            .beg = n - end,
                            .end = n - beg,
                            .src = n - (src + (end - beg))
                        });
                    }
                }
            }

            if (NGV_S[ISA_S[i]] != s) [[likely]] {
                pos_t src = S[SA_S[NGV_S[ISA_S[i]]]];
                pos_t beg = S[i];

                if (!(beg < lst_gr_phr.end &&
                    src - beg == lst_gr_phr.src - lst_gr_phr.beg
                )) {
                    pos_t end = S[i] + LCE_R(S[i], src);

                    if (phr_mode == lpf_lnf_opt && src != 0 && S[i] != 0) {
                        pos_t lce_l = LCE_L(src - 1, S[i] - 1);
                        beg -= lce_l;
                        src -= lce_l;
                    }

                    #ifndef NDEBUG
                    assert(src > beg);

                    for (pos_t j = 0; j < end - beg; j++) {
                        assert(T[src + j] == T[beg + j]);
                    }
                    #endif

                    if (end - beg > 1) [[likely]] {
                        lst_gr_phr = { beg, end, src };
                        
                        LPF[i_p].emplace_back(lpf {
                            .beg = n - end,
                            .end = n - beg,
                            .src = n - (src + (end - beg))
                        });
                    }
                }
            }
        }
    }

    PGV_S.clear();
    NGV_S.clear();
    PGV_S.shrink_to_fit();
    NGV_S.shrink_to_fit();

    if (log) {
        log_phase("lnf", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}