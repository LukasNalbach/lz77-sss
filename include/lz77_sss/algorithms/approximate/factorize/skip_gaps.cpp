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
template <typename output_fnc_t, typename lpf_beg_t, typename next_lpf_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::factorize_skip_gaps(
    output_fnc_t output,
    lpf_beg_t lpf_beg,
    next_lpf_t next_lpf)
{
    if (log) {
        std::cout << "outputting gapped factorization" << std::flush;
    }

    lpf_pos_t lpf_pos = lpf_beg();
    lpf phr;
    lpf phr_nxt = next_lpf(lpf_pos);
    output({ .src = phr_nxt.beg, .len = 0 });

    while (phr_nxt.beg < n) {
        phr = phr_nxt;
        phr_nxt = next_lpf(lpf_pos);
        output({ .src = phr.src, .len = phr.end - phr.beg });

        if (phr_nxt.beg > phr.end) {
            output({ .src = phr_nxt.beg - phr.end, .len = 0 });
        }
    }

    if (log) {
        log_phase("output_gapped", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}