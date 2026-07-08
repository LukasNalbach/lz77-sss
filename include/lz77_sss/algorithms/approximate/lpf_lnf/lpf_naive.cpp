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
void lz77_sss<pos_t>::factorizer<tau, char_t>::build_LPF_naive()
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

        for (uint32_t i = i_min; i < i_max; i++) {
            pos_t src;
            pos_t len = 0;

            if (PSV_S[ISA_S[i]] != s) {
                pos_t src_cur = S[SA_S[PSV_S[ISA_S[i]]]];
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

            if (NSV_S[ISA_S[i]] != s) {
                pos_t src_cur = S[SA_S[NSV_S[ISA_S[i]]]];
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
                LPF[i_p].emplace_back(lpf {
                    .beg = S[i],
                    .end = S[i] + len,
                    .src = src });
            }

            pos_t pos = S[i];

            do {
                i++;
            } while (i < i_max && S[i] < pos + len);
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