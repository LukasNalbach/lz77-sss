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
void lz77_sss<pos_t>::factorizer<tau, char_t>::factorize_greedy(
    output_fnc_t output,
    lpf_beg_t lpf_beg,
    next_lpf_t next_lpf)
{
    if (log) {
        std::cout << "factorizing sequentially" << std::flush;
    }

    lpf_pos_t lpf_pos = lpf_beg();
    lpf p = next_lpf(lpf_pos);

    for (pos_t i = 0; true;) {
        pos_t gap_end = p.beg;

        if (i < gap_end) {
            if (gap_idx.pos() < i) {
                if (i - gap_idx.pos() <= roll_threshold) {
                    do {
                        gap_idx.roll();
                    } while (gap_idx.pos() < i);
                } else {
                    gap_idx.reinit(i);
                }
            }

            do {
                factor f = longest_prev_occ(i);
                i += f.length();

                if (i > gap_end) {
                    if (i <= p.end) {
                        f.len -= i - gap_end;
                        i = gap_end;
                    } else {
                        do {
                            p = next_lpf(lpf_pos);
                        } while (p.end <= i);

                        while (gap_idx.pos() < gap_end) {
                            gap_idx.advance();
                        }

                        gap_end = p.beg;
                    }
                }

                #ifndef NDEBUG
                pos_t i_ = i - f.length();

                assert(f.len == 0 ? f.src == char_to_uchar(T[i_]) : f.src < i_);

                for (pos_t j = 0; j < f.len; j++) {
                    assert(T[i_ + j] == T[f.src + j]);
                }
                #endif

                output(f);
                num_fact++;

                while (gap_idx.pos() < i) {
                    gap_idx.advance();
                }
            } while (i < gap_end);
        }

        if (i == n) {
            break;
        }

        pos_t exc = i - gap_end;

        factor lpf {
            .src = p.src + exc,
            .len = (p.end - p.beg) - exc
        };

        if (gap_idx.pos() == i) {
            factor f = longest_prev_occ(i);

            if (f.len > lpf.len) {
                lpf = f;
            }
        }

        #ifndef NDEBUG
        assert(lpf.src < i);

        for (pos_t j = 0; j < lpf.len; j++) {
            assert(T[i + j] == T[lpf.src + j]);
        }
        #endif

        output(lpf);
        num_fact++;
        i += lpf.len;

        while (p.end <= i) {
            p = next_lpf(lpf_pos);
        }
    }

    if (log) {
        log_phase("factorize", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}