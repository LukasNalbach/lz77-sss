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
void lz77_sss::factorizer<text_t>::factorize_from_lpf(factor_sink& output)
{
    auto lpf_beg = [&]() {
        lpf_cursor_t lpf_cursor { .chunk = 0, .i = 0 };
        while (LPF[lpf_cursor.chunk].empty()) lpf_cursor.chunk++;
        return lpf_cursor;
    };

    auto next_lpf = [&](lpf_cursor_t& lpf_cursor) {
        uint64_t& chunk = lpf_cursor.chunk;
        uint64_t& i = lpf_cursor.i;

        lpf_phrase phr = LPF[chunk][i++];

        if (i == LPF[chunk].size()) [[unlikely]] {
            if (chunk + 1 < LPF.size()) [[likely]] {
                do {
                    chunk++;
                    i = 0;

                    while (i < LPF[chunk].size() && LPF[chunk][i].end <= phr.end) {
                        i++;
                    }

                    if (i < LPF[chunk].size()) {
                        lpf_phrase phr_i = LPF[chunk][i];

                        if (phr_i.beg < phr.end) {
                            uint64_t cut_left = uint64_t(phr.end) - uint64_t(phr_i.beg);
                            phr_i.beg = uint64_t(phr_i.beg) + cut_left;
                            phr_i.src = uint64_t(phr_i.src) + cut_left;
                            LPF[chunk].set(i, phr_i);
                        }
                    }
                } while (i == LPF[chunk].size() && chunk + 1 < LPF.size());
            } else {
                i--;
            }
        }

        return phr;
    };

    if (factorize_gaps_exact) {
        factorize_exact_gaps(output, lpf_beg, next_lpf);
    } else if (fact_mode == auto_gaps) {
        if (gap_ctx > 0) {
            lpf_array trimmed;
            trimmed.init(n + 1);
            lpf_cursor_t lpf_cursor = lpf_beg();
            lpf_phrase phr = next_lpf(lpf_cursor);
            uint64_t prev_end = 0;

            while (phr.beg < n) {
                const lpf_phrase nxt = next_lpf(lpf_cursor);
                const uint64_t beg = phr.beg > prev_end
                    ? std::min<uint64_t>(phr.end, uint64_t(phr.beg) + gap_ctx) : uint64_t(phr.beg);
                const uint64_t end = nxt.beg > phr.end
                    ? std::max<uint64_t>(beg, uint64_t(phr.end) - std::min<uint64_t>(phr.end, gap_ctx)) : uint64_t(phr.end);

                if (end >= beg + 2) {
                    trimmed.emplace_back(lpf_phrase { .beg = beg, .end = end, .src = uint64_t(phr.src) + (beg - uint64_t(phr.beg)) });
                }

                prev_end = phr.end;
                phr = nxt;
            }

            for (lpf_array& phrases : LPF) phrases = lpf_array();
            LPF[0] = std::move(trimmed);
            if (LPF.size() > 1) LPF.back().init(n + 1);
            LPF.back().emplace_back(lpf_phrase { .beg = n, .end = n + 1, .src = 0 });
        }

        factorize_hashed_gaps(output, lpf_beg, next_lpf);
    } else {
        factorize_skip_gaps(output);
    }
}