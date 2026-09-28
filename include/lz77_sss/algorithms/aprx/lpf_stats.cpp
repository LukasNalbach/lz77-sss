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
void lz77_sss::factorizer<text_t>::compute_lpf_stats()
{
    const uint64_t num_chunks = LPF.size();
    std::vector<uint64_t> prev_end(num_chunks + 1, 0);
    uint64_t gaps = 0;
    uint64_t phrases = 0;
    uint64_t len_phrases = 0;

    for (uint64_t c = 0; c < num_chunks; c++) {
        prev_end[c + 1] = LPF[c].empty() ? prev_end[c] : std::max<uint64_t>(prev_end[c], LPF[c].back().end);
    }

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1) reduction(+ : gaps, phrases, len_phrases)
    for (uint64_t c = 0; c < num_chunks; c++) {
        const uint64_t beg = (n * c) / num_chunks;
        const bool continued = c > 0 && prev_end[c] < beg;
        uint64_t b = std::max<uint64_t>(beg, prev_end[c]);
        uint64_t e = (n * (c + 1)) / num_chunks;
        uint64_t i = 0;

        if (!LPF[c].empty()) {
            e = std::max<uint64_t>(e, LPF[c].back().end);
        }

        while (i < LPF[c].size() && LPF[c][i].end <= b) {
            i++;
        }

        phrases += LPF[c].size() - i;

        if (i < LPF[c].size()) {
            len_phrases += LPF[c][i].end - std::max<uint64_t>(LPF[c][i].beg, b);
            if (LPF[c][i].beg > b && !continued) gaps++;
            i++;

            while (i < LPF[c].size()) {
                len_phrases += LPF[c][i].end - LPF[c][i].beg;
                if (LPF[c][i].beg > LPF[c][i - 1].end) gaps++;
                i++;
            }

            if (LPF[c].back().end < e) gaps++;
        } else if (b < e && !continued) {
            gaps++;
        }
    }

    num_gaps += gaps;
    num_lpf_phr += phrases;
    len_lpf_phr += len_phrases;
}