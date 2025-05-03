#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename char_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::build_LPF_opt(std::function<void(uint16_t, lpf&&)> lpf_it)
{
    build_PSV_NSV_S();

    if (log) {
        std::cout << "building LPF" << std::flush;
    }

    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SSA_S = LCE.get_ssa();
    const std::vector<uint32_t>& ISSA_S = LCE.get_issa();
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

        pos_t max_end = b;

        for (uint32_t i = i_min; i < i_max; i++) {
            while (i + 1 < i_max && S[i + 1] <= max_end) {
                i++;
            }

            pos_t lst_end = max_end;
            lpf phr { 0, 0, 0 };

            if (PSV_S[ISSA_S[i]] != s) [[likely]] {
                pos_t src = S[SSA_S[PSV_S[ISSA_S[i]]]];
                pos_t end = S[i] + LCE_R(src, S[i]);

                if (end > lst_end) {
                    pos_t beg = S[i];

                    if (S[i] > lst_end && src != 0 && S[i] != 0) {
                        pos_t lce_l = LCE_L(src - 1, S[i] - 1, S[i] - lst_end);
                        beg -= lce_l;
                        src -= lce_l;
                    }

                    if (beg < lst_end) {
                        pos_t exc = lst_end - beg;
                        beg += exc;
                        src += exc;
                    }

                    if (end > max_end) {
                        max_end = end;
                    }

                    #ifndef NDEBUG
                    assert(src < beg);

                    for (pos_t j = 0; j < end - beg; j++) {
                        assert(T[src + j] == T[beg + j]);
                    }
                    #endif

                    if (end - beg > 1) [[likely]] {
                        phr = { beg, end, src };
                    }
                }
            }

            if (NSV_S[ISSA_S[i]] != s) [[likely]] {
                pos_t src = S[SSA_S[NSV_S[ISSA_S[i]]]];
                pos_t end = S[i] + LCE_R(src, S[i]);

                if (end > lst_end) {
                    pos_t beg = S[i];

                    if (S[i] > lst_end && src != 0 && S[i] != 0) {
                        pos_t lce_l = LCE_L(src - 1, S[i] - 1, S[i] - lst_end);
                        beg -= lce_l;
                        src -= lce_l;
                    }

                    if (beg < lst_end) {
                        pos_t exc = lst_end - beg;
                        beg += exc;
                        src += exc;
                    }

                    if (end > max_end) {
                        max_end = end;
                    }

                    #ifndef NDEBUG
                    assert(src < beg);

                    for (pos_t j = 0; j < end - beg; j++) {
                        assert(T[src + j] == T[beg + j]);
                    }
                    #endif

                    if (end - beg > phr.end - phr.beg) {
                        phr = { beg, end, src };
                    }
                }

                if (phr.end - phr.beg > 1) {
                    lpf_it(i_p, std::move(phr));
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