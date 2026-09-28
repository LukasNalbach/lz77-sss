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

inline lz77_sss::gapped_factorization::gapped_factorization(std::vector<lpf_array>& phrases, uint64_t n)
    : lpf_chunks(&phrases)
    , n(n)
{
    uint64_t prev_end = 0;

    for (uint64_t c = 0; c < phrases.size(); c++) {
        lpf_array& arr = phrases[c];
        uint64_t first = 0;
        while (first < arr.size() && uint64_t(arr[first].end) <= prev_end) first++;
        if (first == arr.size()) continue;
        lpf_phrase phr = arr[first];

        if (uint64_t(phr.beg) < prev_end) {
            const uint64_t offs = prev_end - uint64_t(phr.beg);
            phr.beg = prev_end;
            phr.src = uint64_t(phr.src) + offs;
            arr.set(first, phr);
        }

        sections.emplace_back(sect_t { .chunk = c, .first = first, .beg = uint64_t(phr.beg), .nxt = n });
        prev_end = std::max<uint64_t>(prev_end, uint64_t(arr.back().end));
    }

    for (uint64_t s = 0; s + 1 < sections.size(); s++) sections[s].nxt = sections[s + 1].beg;
}

inline void lz77_sss::gapped_factorization::emit_section(uint64_t s, factor_sink& out) const
{
    const sect_t& sc = sections[s];
    const lpf_array& arr = (*lpf_chunks)[sc.chunk];
    if (s == 0) out(factor::gap(sc.beg));

    for (uint64_t k = sc.first; k < arr.size(); k++) {
        const lpf_phrase phr = arr[k];
        if (uint64_t(phr.beg) >= n) break;
        out(factor { .src = phr.src, .len = uint64_t(phr.end) - uint64_t(phr.beg) });
        const uint64_t nxt = k + 1 < arr.size() ? uint64_t(arr[k + 1].beg) : sc.nxt;
        if (nxt > uint64_t(phr.end)) out(factor::gap(nxt - uint64_t(phr.end)));
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::factorize_skip_gaps(factor_sink& output)
{
    if (log) {
        std::cout << "outputting gapped factorization" << std::flush;
    }

    gapped_factorization gapped(LPF, n);

    if (gapped_fnc) {
        gapped_fnc(gapped);
    } else {
        for (uint64_t s = 0; s < gapped.num_sections(); s++) gapped.emit_section(s, output);
    }

    if (log) {
        record_phase_time("output_gapped", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}
