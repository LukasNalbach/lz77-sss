#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::build_LPF_greedy() {
    build_PSV_NSV_S();

    if (log) {
        std::cout << "building LPF" << std::flush;
    }
    
    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SSA_S = LCE.get_ssa();
    const std::vector<uint32_t>& ISSA_S = LCE.get_issa();
    pos_t s = S.size();

    for (uint32_t i = 0; i < s;) {
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
                src = src_cur;
                len = len_cur;
            }
        }

        if (len > 0) {
            LPF.emplace_back(lpf {
                .beg = S[i],
                .end = S[i] + len,
                .src = src
            });
        }

        pos_t pos = S[i];

        do {
            i++;
        } while (i < s && S[i] < pos + len);
    }

    NSV_S.clear();
    PSV_S.clear();
    NSV_S.shrink_to_fit();
    PSV_S.shrink_to_fit();

    if (log) {
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::build_LNF_greedy() {
    build_PGV_NGV_S();
    
    if (log) {
        std::cout << "building LNF" << std::flush;
    }
    
    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SSA_S = LCE.get_ssa();
    const std::vector<uint32_t>& ISSA_S = LCE.get_issa();
    pos_t s = S.size();

    for (uint32_t i = 0; i < s;) {
        pos_t src;
        pos_t len = 0;

        if (PGV_S[ISSA_S[i]] > 0) {
            pos_t src_cur = S[SSA_S[PGV_S[ISSA_S[i]]]];
            pos_t len_cur = LCE_R(src_cur, S[i]);

            #ifndef NDEBUG
            assert(src_cur > S[i]);

            for (pos_t j = 0; j < len_cur; j++) {
                assert(T[src_cur + j] == T[S[i] + j]);
            }
            #endif

            if (len_cur > len) {
                src = src_cur;
                len = len_cur;
            }
        }

        if (NGV_S[ISSA_S[i]] < s) {
            pos_t src_cur = S[SSA_S[NGV_S[ISSA_S[i]]]];
            pos_t len_cur = LCE_R(src_cur, S[i]);

            #ifndef NDEBUG
            assert(src_cur > S[i]);
            
            for (pos_t j = 0; j < len_cur; j++) {
                assert(T[src_cur + j] == T[S[i] + j]);
            }
            #endif

            if (len_cur > len) {
                src = src_cur;
                len = len_cur;
            }
        }

        if (len > 0) {
            LPF.emplace_back(lpf {
                .beg = n - S[i] - len,
                .end = n - S[i],
                .src = n - src - len
            });
        }
        
        pos_t pos = S[i];

        do {
            i++;
        } while (i < s && S[i] < pos + len);
    }

    PGV_S.clear();
    NGV_S.clear();
    PGV_S.shrink_to_fit();
    NGV_S.shrink_to_fit();

    if (log) {
        time = log_runtime(time);
    }
}