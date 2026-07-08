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

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename char_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::build_LPF_opt()
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

        pos_t max_end = b;

        for (uint32_t i = i_min; i < i_max; i++) {
            while (i + 1 < i_max && S[i + 1] <= max_end) {
                i++;
            }

            pos_t lst_end = max_end;
            lpf phr { 0, 0, 0 };

            if (PSV_S[ISA_S[i]] != s) [[likely]] {
                pos_t src = S[SA_S[PSV_S[ISA_S[i]]]];
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

            if (NSV_S[ISA_S[i]] != s) [[likely]] {
                pos_t src = S[SA_S[NSV_S[ISA_S[i]]]];
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
                    LPF[i_p].emplace_back(phr);
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