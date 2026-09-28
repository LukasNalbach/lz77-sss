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

template <typename text_t>
void lz77_sss::factorizer<text_t>::build_LPF()
{
    if (log) {
        std::cout << "building LPF" << std::flush;
    }

    const auto& S = LCE.get_sync_set();
    const auto& SA_S = LCE.get_sa_s();
    const auto& ISA_S = LCE.get_isa_s();
    uint64_t s = S.size();
    const min_tree<bit_aligned_vector> smaller(SA_S, s, p);
    const uint64_t num_chunks = std::clamp<uint64_t>(s, 1, uint64_t(p) * lpf_chunks_per_thr);
    LPF.clear();
    LPF.resize(num_chunks);
    for (lpf_array& phrases : LPF) phrases.init(n + 1);

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t chunk = 0; chunk < num_chunks; chunk++) {
        uint64_t b = (n * chunk) / num_chunks;
        uint64_t e = (n * (chunk + 1)) / num_chunks;

        uint32_t i_min = bin_search_min_geq<uint64_t, uint32_t>(
            true, 0, s, [&](uint32_t i) { return i == s || S[i] >= b; });
        uint32_t i_max = bin_search_min_geq<uint64_t, uint32_t>(
            true, 0, s, [&](uint32_t i) { return i == s || S[i] >= e; });

        uint64_t max_end = b;

        for (uint32_t i = i_min; i < i_max; i++) {
            while (i + 1 < i_max && S[i + 1] <= max_end) {
                i++;
            }

            uint64_t prev_end = max_end;
            lpf_phrase phr { 0, 0, 0 };

            const uint64_t x = ISA_S[i];
            const uint64_t psv = smaller.previous_smaller(x);

            if (psv != s) [[likely]] {
                uint64_t src = S[SA_S[psv]];
                uint64_t end = S[i] + LCE_R(src, S[i]);

                if (end > prev_end) {
                    uint64_t beg = S[i];

                    if (S[i] > prev_end && src != 0 && S[i] != 0) {
                        uint64_t lce_l = LCE_L(src - 1, S[i] - 1, S[i] - prev_end);
                        beg -= lce_l;
                        src -= lce_l;
                    }

                    if (beg < prev_end) {
                        uint64_t cut_left = prev_end - beg;
                        beg += cut_left;
                        src += cut_left;
                    }

                    if (end > max_end) {
                        max_end = end;
                    }

                    #ifndef NDEBUG
                    assert(src < beg);

                    for (uint64_t j = 0; j < end - beg; j++) {
                        assert(T[src + j] == T[beg + j]);
                    }
                    #endif

                    if (end - beg > 1) [[likely]] {
                        phr = { beg, end, src };
                    }
                }
            }

            const uint64_t nsv = smaller.next_smaller(x);

            if (nsv != s) [[likely]] {
                uint64_t src = S[SA_S[nsv]];
                uint64_t end = S[i] + LCE_R(src, S[i]);

                if (end > prev_end) {
                    uint64_t beg = S[i];

                    if (S[i] > prev_end && src != 0 && S[i] != 0) {
                        uint64_t lce_l = LCE_L(src - 1, S[i] - 1, S[i] - prev_end);
                        beg -= lce_l;
                        src -= lce_l;
                    }

                    if (beg < prev_end) {
                        uint64_t cut_left = prev_end - beg;
                        beg += cut_left;
                        src += cut_left;
                    }

                    if (end > max_end) {
                        max_end = end;
                    }

                    #ifndef NDEBUG
                    assert(src < beg);

                    for (uint64_t j = 0; j < end - beg; j++) {
                        assert(T[src + j] == T[beg + j]);
                    }
                    #endif

                    if (end - beg > phr.end - phr.beg) {
                        phr = { beg, end, src };
                    }
                }
            }

            if (phr.end - phr.beg > 1) {
                LPF[chunk].emplace_back(phr);
            }
        }
    }

    if (log) {
        record_phase_time("lpf", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}