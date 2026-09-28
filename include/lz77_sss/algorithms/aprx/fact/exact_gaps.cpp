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
template <typename lpf_beg_t, typename next_lpf_t>
void lz77_sss::factorizer<text_t>::factorize_exact_gaps(
    factor_sink& output,
    lpf_beg_t lpf_beg,
    next_lpf_t next_lpf)
{
    if (log) {
        std::cout << "concatenating gaps" << std::flush;
    }

    std::vector<uint64_t> reg_beg;
    std::vector<uint64_t> reg_end;
    std::vector<lpf_phrase> phrases;

    {
        lpf_cursor_t lpf_cursor = lpf_beg();
        uint64_t prev_end = 0;

        while (true) {
            const lpf_phrase phr = next_lpf(lpf_cursor);
            const uint64_t beg = std::min<uint64_t>(phr.beg, n);

            if (beg > prev_end) {
                const uint64_t ctx_beg = prev_end > gap_ctx ? prev_end - gap_ctx : 0;
                const uint64_t ctx_end = std::min<uint64_t>(n, beg + gap_ctx);

                if (!reg_end.empty() && ctx_beg <= reg_end.back()) {
                    reg_end.back() = std::max<uint64_t>(reg_end.back(), ctx_end);
                } else {
                    reg_beg.emplace_back(ctx_beg);
                    reg_end.emplace_back(ctx_end);
                }
            }

            if (phr.beg >= n) break;
            if (gap_ctx > 0) phrases.emplace_back(phr);
            prev_end = phr.end;
        }
    }

    const uint64_t num_regs = reg_beg.size();
    std::vector<uint64_t> reg_off(num_regs + 1, 0);

    for (uint64_t r = 0; r < num_regs; r++) {
        reg_off[r + 1] = reg_off[r] + (reg_end[r] - reg_beg[r]) + 1;
    }

    const uint64_t len_G = reg_off[num_regs];
    std::string G;
    no_init_resize(G, len_G);
    std::array<uint64_t, 256> hist { };

    #pragma omp parallel num_threads(p)
    {
        std::array<uint64_t, 256> hist_thr { };

        #pragma omp for schedule(dynamic, 256)
        for (uint64_t r = 0; r < num_regs; r++) {
            auto cursor = T.cursor_at(reg_beg[r]);
            char* dst = G.data() + reg_off[r];

            for (uint64_t i = 0; i < reg_end[r] - reg_beg[r]; i++) {
                const uint8_t c = cursor.next();
                dst[i] = char(c);
                hist_thr[c]++;
            }
        }

        #pragma omp critical
        for (uint64_t c = 0; c < 256; c++) hist[c] += hist_thr[c];
    }

    uint8_t sep = 0;

    for (uint64_t c = 256; c-- > 0;) {
        if (hist[c] == 0) {
            sep = uint8_t(c);
            break;
        }
    }

    #pragma omp parallel for num_threads(p)
    for (uint64_t r = 0; r < num_regs; r++) {
        G[reg_off[r + 1] - 1] = char(sep);
    }

    if (log) {
        std::cout << " (" << num_regs << " gaps, " << format_size(len_G) << ")";
        time = log_runtime(time);
        std::cout << "factorizing gaps (exact)" << std::flush;
    }

    const direct_text G_text(G.data(), len_G);
    const uint64_t num_blks = std::min<uint64_t>(num_regs, uint64_t { p } * par_sects_per_thr);
    std::vector<std::vector<factor>> facts(num_blks);
    std::vector<uint64_t> num_facts(num_regs);

    auto factorize_gaps = [&](const auto& idx) {
        #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
        for (uint64_t b = 0; b < num_blks; b++) {
            const uint64_t blk_r_beg = bin_search_min_geq<uint64_t, uint64_t>(
                (len_G * b) / num_blks, 0, num_regs, [&](uint64_t r) { return reg_off[r]; });
            const uint64_t blk_r_end = bin_search_min_geq<uint64_t, uint64_t>(
                (len_G * (b + 1)) / num_blks, 0, num_regs, [&](uint64_t r) { return reg_off[r]; });
            uint64_t j = gap_ctx == 0 || blk_r_beg == blk_r_end ? 0 : bin_search_min_geq<uint64_t, uint64_t>(
                reg_beg[blk_r_beg] + 1, 0, phrases.size(), [&](uint64_t i) { return i == phrases.size() ? n + 1 : uint64_t(phrases[i].end); });

            for (uint64_t r = blk_r_beg; r < blk_r_end; r++) {
                const uint64_t end_G = reg_off[r + 1] - 1;
                const uint64_t facts_beg = facts[b].size();

                for (uint64_t x = reg_off[r]; x < end_G;) {
                    const uint64_t pos = reg_beg[r] + (x - reg_off[r]);
                    const auto occ = idx.lpf(x, end_G - x);
                    factor f = factor::literal(T.to_char(T[pos]));

                    if (occ.len > 0) {
                        const uint64_t r_src = bin_search_max_leq<uint64_t, uint64_t>(
                            occ.src, 0, num_regs - 1, [&](uint64_t i) { return reg_off[i]; });
                        const uint64_t src = reg_beg[r_src] + (occ.src - reg_off[r_src]);
                        const uint64_t len = LCE_R(src, pos);
                        if (len > 0) f = factor { .src = src, .len = len };
                    }

                    if (gap_ctx > 0) {
                        while (j < phrases.size() && phrases[j].end <= pos) j++;

                        if (j < phrases.size() && phrases[j].beg <= pos) {
                            const uint64_t src = uint64_t(phrases[j].src) + (pos - uint64_t(phrases[j].beg));
                            const uint64_t len = LCE_R(src, pos);
                            if (len > f.len || (len == f.len && src > f.src)) f = factor { .src = src, .len = len };
                        }
                    }

                    facts[b].emplace_back(f);
                    x += f.text_len();
                }

                num_facts[r] = facts[b].size() - facts_beg;
            }
        }
    };

    if (len_G <= INT32_MAX) {
        factorize_gaps(exact_gap_index<direct_text, int32_t>(G_text, reinterpret_cast<const uint8_t*>(G.data()), len_G, p));
    } else {
        factorize_gaps(exact_gap_index<direct_text, sa_int40_t>(G_text, reinterpret_cast<const uint8_t*>(G.data()), len_G, p));
    }

    G = std::string();
    phrases = std::vector<lpf_phrase>();

    if (log) {
        time = log_runtime(time);
        std::cout << "interleaving factors" << std::flush;
    }

    uint64_t out_end = 0;
    merging_output merged(output, num_fact);

    auto emit = [&](uint64_t pos, factor f) {
        uint64_t len = f.text_len();
        if (pos + len <= out_end) return;

        if (pos < out_end) {
            const uint64_t cut_left = out_end - pos;
            f.src = f.src + cut_left;
            f.len = f.len - cut_left;
            len -= cut_left;
            pos = out_end;
        }

        #ifndef NDEBUG
        assert((f.is_literal() && f.src == T.to_char(T[pos])) || f.src < pos);
        for (uint64_t x = 0; x < f.len; x++) assert(T[f.src + x] == T[pos + x]);
        #endif

        merged.add(pos, f);
        out_end = pos + len;
    };

    uint64_t b = 0;
    uint64_t i = 0;

    auto emit_region = [&](uint64_t r) {
        uint64_t pos = reg_beg[r];

        for (uint64_t c = 0; c < num_facts[r]; c++) {
            while (i == facts[b].size()) {
                facts[b] = std::vector<factor>();
                b++;
                i = 0;
            }

            const factor f = facts[b][i++];
            emit(pos, f);
            pos += f.text_len();
        }
    };

    lpf_cursor_t lpf_cursor = lpf_beg();
    uint64_t r = 0;

    while (true) {
        const lpf_phrase phr = next_lpf(lpf_cursor);
        if (phr.beg >= n) break;
        uint64_t beg = phr.beg;

        while (r < num_regs && reg_beg[r] <= beg) {
            emit_region(r);
            beg = std::max<uint64_t>(beg, reg_end[r]);
            r++;
        }

        const uint64_t end = r < num_regs ? std::min<uint64_t>(phr.end, reg_beg[r]) : uint64_t(phr.end);
        if (end > beg) emit(beg, factor { .src = uint64_t(phr.src) + (beg - uint64_t(phr.beg)), .len = end - beg });
    }

    while (r < num_regs) emit_region(r++);
    merged.flush();

    if (log) {
        record_phase_time("factorize", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}
