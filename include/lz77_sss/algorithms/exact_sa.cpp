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
template <typename sym_t>
std::vector<sym_t> lz77_sss::factorizer<text_t>::symbols_of(const text_t& text, uint16_t p)
{
    const uint64_t size = text.size();
    std::vector<sym_t> symbols;
    no_init_resize(symbols, size);

    parallel_chunks(size, p, [&](uint64_t beg, uint64_t end) {
        auto cursor = text.cursor_at(beg);
        for (uint64_t i = beg; i < end; i++) symbols[i] = sym_t(cursor.next());
    });

    return symbols;
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::factorize_exact_sa(factor_sink& output)
{
    if (n <= INT32_MAX) {
        factorize_exact_sa<int32_t>(output);
    } else {
        factorize_exact_sa<sa_int40_t>(output);
    }
}

template <typename text_t>
template <typename sa_t>
void lz77_sss::factorizer<text_t>::factorize_exact_sa(factor_sink& output)
{
    if (log) {
        std::cout << "building SA" << std::flush;
    }

    std::optional<exact_gap_index<text_t, sa_t>> idx;
    bool sa_built = false;

    auto step = [&](uint64_t bytes) {
        if (!log) return;
        record_phase_time(sa_built ? "isa" : "sa", time_diff_ns(time, now()));
        std::cout << " (" << format_size(bytes) << ")";
        time = log_runtime(time);
        std::cout << (sa_built ? "computing the exact factorization" : "building ISA") << std::flush;
        sa_built = true;
    };

    if constexpr (std::is_same_v<text_t, direct_text>) {
        idx.emplace(T, reinterpret_cast<const uint8_t*>(T.data()), n, p, step);
    } else if constexpr (text_t::is_byte_text) {
        idx.emplace(T, symbols_of<uint8_t>(T, p), n, p, step);
    } else if constexpr (std::is_same_v<sa_t, int32_t>) {
        if constexpr (std::is_same_v<text_t, int_direct_text>) {
            idx.emplace(T, reinterpret_cast<int32_t*>(T.data()), n, T.sigma(), p, step);
        } else {
            idx.emplace(T, symbols_of<int32_t>(T, p), n, T.sigma(), p, step);
        }
    } else {
        idx.emplace(T, symbols_of<sa_int40_t>(T, p), n, T.sigma(), p, step);
    }

    auto factor_at = [&](uint64_t i) {
        const auto occ = idx->lpf(i, n - i);
        if (occ.len == 0) return factor::literal(T.to_char(T[i]));
        return factor { .src = occ.src, .len = occ.len };
    };

    const uint64_t num_sect = p == 1 ? 1 : std::clamp<uint64_t>(n / min_sa_sect_len, 1,
        std::max<uint64_t>(uint64_t { p } * par_sects_per_thr, n / max_sa_sect_len));
    auto sect_beg = [&](uint64_t k) { return (n * k) / num_sect; };
    std::vector<std::vector<factor>> facts(num_sect);
    num_fact = 0;
    uint64_t pos = 0;

    auto emit = [&](factor f) {
        #ifndef NDEBUG
        assert((f.is_literal() && f.src == T.to_char(T[pos])) || f.src < pos);
        for (uint64_t x = 0; x < f.len; x++) assert(T[f.src + x] == T[pos + x]);
        #endif

        output(f);
        num_fact++;
        pos += f.text_len();
    };

    auto factorize_sect = [&](uint64_t k) {
        const uint64_t e = sect_beg(k + 1);

        for (uint64_t i = sect_beg(k); i < e;) {
            const factor f = factor_at(i);
            facts[k].emplace_back(f);
            i += f.text_len();
        }
    };

    auto emit_sect = [&](uint64_t k) {
        uint64_t q = sect_beg(k);

        for (const factor& f : facts[k]) {
            const uint64_t beg = q;
            q += f.text_len();
            while (pos < beg) emit(factor_at(pos));
            if (beg == pos) emit(f);
        }

        facts[k] = std::vector<factor>();
    };

    if (p > 1) {
        std::vector<std::atomic<uint8_t>> done(num_sect);
        std::atomic<uint64_t> next_sect = 0;
        std::atomic<uint64_t> sects_emitted = 0;
        const uint64_t max_ahead = uint64_t { p } * sa_sects_ahead_per_thr;

        #pragma omp parallel num_threads(p)
        {
            if (omp_get_thread_num() == 0) {
                const bool alone = omp_get_num_threads() == 1;

                for (uint64_t k = 0; k < num_sect; k++) {
                    if (alone) factorize_sect(k);
                    spin_until([&]() { return alone || done[k].load(std::memory_order_acquire) != 0; });
                    emit_sect(k);
                    sects_emitted.store(k + 1, std::memory_order_release);
                }
            } else {
                for (uint64_t k; (k = next_sect.fetch_add(1, std::memory_order_relaxed)) < num_sect;) {
                    spin_until([&]() { return k < sects_emitted.load(std::memory_order_acquire) + max_ahead; });
                    factorize_sect(k);
                    done[k].store(1, std::memory_order_release);
                }
            }
        }
    }

    while (pos < n) emit(factor_at(pos));

    if (log) {
        record_phase_time("factorize_sa", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}
