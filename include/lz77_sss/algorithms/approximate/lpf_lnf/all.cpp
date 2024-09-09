#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::build_LPF_all(std::function<void(lpf&&, pos_t)> lpf_it) {
    build_PSV_NSV_S();

    if (log) {
        std::cout << "building LPF" << std::flush;
    }

    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SSA_S = LCE.get_ssa();
    const std::vector<uint32_t>& ISSA_S = LCE.get_issa();
    pos_t s = S.size();

    for (uint32_t i = 0; i < s; i++) {
        if (PSV_S[ISSA_S[i]] > 0) [[likely]] {
            pos_t src = S[SSA_S[PSV_S[ISSA_S[i]]]];
            pos_t lce_l = 0;

            if (src != 0 && S[i] != 0) [[likely]] {
                lce_l = LCE_L(src - 1, S[i] - 1, tau);
            }
            
            pos_t end = S[i] + LCE_R(src, S[i]);
            pos_t beg = S[i] - lce_l;
            src -= lce_l;

            #ifndef NDEBUG
            assert(src < beg);

            for (pos_t j = 0; j < end - beg; j++) {
                assert(T[src + j] == T[beg + j]);
            }
            #endif

            if (end - beg > 1) [[likely]] {
                lpf_it(lpf {beg, end, src}, S[i]);
            }
        }

        if (NSV_S[ISSA_S[i]] < s) [[likely]] {
            pos_t src = S[SSA_S[NSV_S[ISSA_S[i]]]];
            pos_t lce_l = 0;

            if (src != 0 && S[i] != 0) [[likely]] {
                lce_l = LCE_L(src - 1, S[i] - 1, tau);
            }

            pos_t end = S[i] + LCE_R(src, S[i]);
            pos_t beg = S[i] - lce_l;
            src -= lce_l;

            #ifndef NDEBUG
            assert(src < beg);
            
            for (pos_t j = 0; j < end - beg; j++) {
                assert(T[src + j] == T[beg + j]);
            }
            #endif

            if (end - beg > 1) [[likely]] {
                lpf_it(lpf {beg, end, src}, S[i]);
            }
        }
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
void lz77_sss<pos_t>::factorizer<tau>::build_LNF_all(std::function<void(lpf&&, pos_t)> lpf_it) {
    build_PGV_NGV_S();

    if (log) {
        std::cout << "building LNF" << std::flush;
    }
    
    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SSA_S = LCE.get_ssa();
    const std::vector<uint32_t>& ISSA_S = LCE.get_issa();
    pos_t s = S.size();

    for (uint32_t i = 0; i < s; i++) {
        if (PGV_S[ISSA_S[i]] > 0) [[likely]] {
            pos_t src = S[SSA_S[PGV_S[ISSA_S[i]]]];
            pos_t lce_l = 0;

            if (src != 0 && S[i] != 0) [[likely]] {
                lce_l = LCE_L(src - 1, S[i] - 1, tau);
            }

            pos_t end = S[i] + LCE_R(src, S[i]);
            pos_t beg = S[i] - lce_l;
            src -= lce_l;

            #ifndef NDEBUG
            assert(src > beg);

            for (pos_t j = 0; j < end - beg; j++) {
                assert(T[src + j] == T[beg + j]);
            }
            #endif

            if (end - beg > 1) [[likely]] {
                lpf_it(lpf {
                    .beg = n - end,
                    .end = n - beg,
                    .src = n - (src + (end - beg))
                }, S[i]);
            }
        }

        if (NGV_S[ISSA_S[i]] < s) [[likely]] {
            pos_t src = S[SSA_S[NGV_S[ISSA_S[i]]]];
            pos_t lce_l = 0;

            if (src != 0 && S[i] != 0) [[likely]] {
                lce_l = LCE_L(src - 1, S[i] - 1, tau);
            }

            pos_t end = S[i] + LCE_R(src, S[i]);
            pos_t beg = S[i] - lce_l;
            src -= lce_l;

            #ifndef NDEBUG
            assert(src > beg);
            
            for (pos_t j = 0; j < end - beg; j++) {
                assert(T[src + j] == T[beg + j]);
            }
            #endif

            if (end - beg > 1) [[likely]] {
                lpf_it(lpf {
                    .beg = n - end,
                    .end = n - beg,
                    .src = n - (src + (end - beg))
                }, S[i]);
            }
        }
    }

    PGV_S.clear();
    NGV_S.clear();
    PGV_S.shrink_to_fit();
    NGV_S.shrink_to_fit();

    if (log) {
        time = log_runtime(time);
    }
}