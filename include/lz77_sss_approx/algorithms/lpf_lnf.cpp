#pragma once

#include <lz77_sss_approx/lz77_sss_approx.hpp>

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
template <lz77_sss_approx<pos_t>::lpf_mode mode>
void lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::build_LPF_S() {
    if (log) {
        std::cout << "building LPF" << std::flush;
    }
    
    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SSA_S = LCE.get_ssa();
    const std::vector<uint32_t>& ISSA_S = LCE.get_issa();
    pos_t s = S.size();

    #ifndef NDEBUG
    for (uint32_t i = 1; i < s; i++) {
        assert(S[i - 1] < S[i]);
    }
    #endif

    std::vector<uint32_t> PSV_S;
    std::vector<uint32_t> NSV_S;
    no_init_resize(PSV_S, s);
    no_init_resize(NSV_S, s);

    PSV_S[0] = 0;
    NSV_S[s - 1] = s;

    for (uint32_t i = 1; i < s; i++) {
        uint32_t j = i - 1;

        while (SSA_S[j] > SSA_S[i] && j > 0) {
            j = PSV_S[j];
        }

        PSV_S[i] = j;
        
        #ifndef NDEBUG
        assert(PSV_S[i] == 0 || S[SSA_S[PSV_S[i]]] < S[SSA_S[i]]);
        #endif
    }

    for (int64_t i = s - 1; i >= 0; i--) {
        uint32_t j = i + 1;

        while (SSA_S[j] > SSA_S[i] && j < s) {
            j = NSV_S[j];
        }

        NSV_S[i] = j;
        
        #ifndef NDEBUG
        assert(NSV_S[i] == s || S[SSA_S[NSV_S[i]]] < S[SSA_S[i]]);
        #endif
    }

    if constexpr (mode == naive) {
        LPF.reserve(s);

        for (uint32_t i = 0; i < s;) {
            pos_t pos;
            pos_t src;
            pos_t len = 0;

            if (PSV_S[ISSA_S[i]] > 0) {
                pos_t src_cur = S[SSA_S[PSV_S[ISSA_S[i]]]];
                pos_t len_cur = LCE_R(src_cur, S[i]);

                #ifndef NDEBUG
                assert(src_cur < S[i]);

                for (pos_t j = 0; j < len_cur; j++) {
                    assert(T[src_cur + j] == T[S[i] + j]);
                }
                #endif

                if (len_cur > len) {
                    pos = S[i];
                    src = src_cur;
                    len = len_cur;
                }
            }

            if (NSV_S[ISSA_S[i]] < s) {
                pos_t src_cur = S[SSA_S[NSV_S[ISSA_S[i]]]];
                pos_t len_cur = LCE_R(src_cur, S[i]);

                #ifndef NDEBUG
                assert(src_cur < S[i]);
                
                for (pos_t j = 0; j < len_cur; j++) {
                    assert(T[src_cur + j] == T[S[i] + j]);
                }
                #endif

                if (len_cur > len) {
                    pos = S[i];
                    src = src_cur;
                    len = len_cur;
                }
            }

            if (len > 0) {
                LPF.emplace_back(lpf {
                    .beg = pos,
                    .end = pos + len,
                    .src = src
                });
            }

            do {
                i++;
            } while (i < s && S[i] < pos + len);
        }
    } else if (mode == all_extend_left) {
        LPF.reserve(2 * s);

        for (uint32_t i = 0; i < s; i++) {
            if (PSV_S[ISSA_S[i]] > 0) {
                pos_t src = S[SSA_S[PSV_S[ISSA_S[i]]]];
                pos_t lce_l = src == 0 ? 0 : LCE_L(src - 1, S[i] - 1);
                pos_t end = S[i] + LCE_R(src, S[i]);
                pos_t beg = S[i] - lce_l;
                src -= lce_l;

                #ifndef NDEBUG
                assert(src < beg);

                for (pos_t j = 0; j < end - beg; j++) {
                    assert(T[src + j] == T[beg + j]);
                }
                #endif

                if (end - beg > 1) {
                    LPF.emplace_back(lpf {beg, end, src});
                }
            }

            if (NSV_S[ISSA_S[i]] < s) {
                pos_t src = S[SSA_S[NSV_S[ISSA_S[i]]]];
                pos_t lce_l = src == 0 ? 0 : LCE_L(src - 1, S[i] - 1);
                pos_t end = S[i] + LCE_R(src, S[i]);
                pos_t beg = S[i] - lce_l;
                src -= lce_l;

                #ifndef NDEBUG
                assert(src < beg);
                
                for (pos_t j = 0; j < end - beg; j++) {
                    assert(T[src + j] == T[beg + j]);
                }
                #endif

                if (end - beg > 1) {
                    LPF.emplace_back(lpf {beg, end, src});
                }
            }
        }
    }

    if (log) {
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
template <lz77_sss_approx<pos_t>::lpf_mode mode>
void lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::build_LNF_S() {
    if (log) {
        std::cout << "building LNF" << std::flush;
    }
    
    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SSA_S = LCE.get_ssa();
    const std::vector<uint32_t>& ISSA_S = LCE.get_issa();
    pos_t s = S.size();

    #ifndef NDEBUG
    for (uint32_t i = 1; i < s; i++) {
        assert(S[i-1] < S[i]);
    }
    #endif

    std::vector<uint32_t> PGV_S;
    std::vector<uint32_t> NGV_S;
    no_init_resize(PGV_S, s);
    no_init_resize(NGV_S, s);

    PGV_S[0] = 0;
    NGV_S[s - 1] = s;

    for (uint32_t i = 1; i < s; i++) {
        uint32_t j = i - 1;

        while (SSA_S[j] < SSA_S[i] && j > 0) {
            j = PGV_S[j];
        }

        PGV_S[i] = j;
        
        #ifndef NDEBUG
        assert(PGV_S[i] == 0 || S[SSA_S[PGV_S[i]]] > S[SSA_S[i]]);
        #endif
    }

    for (int64_t i = s - 1; i >= 0; i--) {
        uint32_t j = i + 1;

        while (SSA_S[j] < SSA_S[i] && j < s) {
            j = NGV_S[j];
        }

        NGV_S[i] = j;
        
        #ifndef NDEBUG
        assert(NGV_S[i] == s || S[SSA_S[NGV_S[i]]] > S[SSA_S[i]]);
        #endif
    }

    if constexpr (mode == naive) {
        LPF.reserve(s);

        for (uint32_t i = 0; i < s;) {
            pos_t pos;
            pos_t src;
            pos_t len = 0;

            if (PGV_S[ISSA_S[i]] > 0) {
                pos_t src_cur = S[SSA_S[PGV_S[ISSA_S[i]]]];
                pos_t lce_l = src_cur == 0 ? 0 : LCE_L(src_cur - 1, S[i] - 1);
                pos_t len_cur = LCE_R(src_cur, S[i]) + lce_l;
                pos_t pos_cur = S[i] - lce_l;
                src_cur -= lce_l;

                #ifndef NDEBUG
                assert(src_cur > pos_cur);

                for (pos_t j = 0; j < len_cur; j++) {
                    assert(T[src_cur + j] == T[pos_cur + j]);
                }
                #endif

                if (len_cur > len) {
                    pos = pos_cur;
                    src = src_cur;
                    len = len_cur;
                }
            }

            if (NGV_S[ISSA_S[i]] < s) {
                pos_t src_cur = S[SSA_S[NGV_S[ISSA_S[i]]]];
                pos_t lce_l = src_cur == 0 ? 0 : LCE_L(src_cur - 1, S[i] - 1);
                pos_t len_cur = LCE_R(src_cur, S[i]) + lce_l;
                pos_t pos_cur = S[i] - lce_l;
                src_cur -= lce_l;

                #ifndef NDEBUG
                assert(src_cur > pos_cur);
                
                for (pos_t j = 0; j < len_cur; j++) {
                    assert(T[src_cur + j] == T[pos_cur + j]);
                }
                #endif

                if (len_cur > len) {
                    pos = pos_cur;
                    src = src_cur;
                    len = len_cur;
                }
            }

            if (len > 0) {
                LPF.emplace_back(lpf {
                    .beg = n - (pos + len),
                    .end = n - pos,
                    .src = n - (src + len)
                });
            }
            
            do {
                i++;
            } while (i < s && S[i] < pos + len);
        }
    } else if (mode == all_extend_left) {
        LPF.reserve(2 * s);

        for (uint32_t i = 0; i < s; i++) {
            if (PGV_S[ISSA_S[i]] > 0) {
                pos_t src = S[SSA_S[PGV_S[ISSA_S[i]]]];
                pos_t lce_l = src == 0 ? 0 : LCE_L(src - 1, S[i] - 1);
                pos_t end = S[i] + LCE_R(src, S[i]);
                pos_t beg = S[i] - lce_l;
                src -= lce_l;

                #ifndef NDEBUG
                assert(src < beg);

                for (pos_t j = 0; j < end - beg; j++) {
                    assert(T[src + j] == T[beg + j]);
                }
                #endif

                if (end - beg > 1) {
                    LPF.emplace_back(lpf {
                        .beg = n - end,
                        .end = n - beg,
                        .src = n - (src + (end - beg))
                    });
                }
            }

            if (NGV_S[ISSA_S[i]] < s) {
                pos_t src = S[SSA_S[NGV_S[ISSA_S[i]]]];
                pos_t lce_l = src == 0 ? 0 : LCE_L(src - 1, S[i] - 1);
                pos_t end = S[i] + LCE_R(src, S[i]);
                pos_t beg = S[i] - lce_l;
                src -= lce_l;

                #ifndef NDEBUG
                assert(src < beg);
                
                for (pos_t j = 0; j < end - beg; j++) {
                    assert(T[src + j] == T[beg + j]);
                }
                #endif

                if (end - beg > 1) {
                    LPF.emplace_back(lpf {
                        .beg = n - end,
                        .end = n - beg,
                        .src = n - (src + (end - beg))
                    });
                }
            }
        }
    }

    if (log) {
        time = log_runtime(time);
    }
}