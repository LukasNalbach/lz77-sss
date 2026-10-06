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

template <std::input_iterator fact_it_t, typename out_it_t>
void lz77_sss::decode(fact_it_t fact_it, out_it_t out_it, uint64_t output_size)
{
    using char_t = std::iter_value_t<out_it_t>;
    factor f;
    uint64_t pos = 0;

    while (pos < output_size) {
        f = *fact_it++;

        if (f.is_literal()) {
            if constexpr (sizeof(char_t) == 1) {
                out_it[pos] = unsigned_to_char<uint64_t, char_t>(f.src);
            } else {
                out_it[pos] = char_t(uint64_t(f.src));
            }

            pos++;
        } else {
            #ifndef NDEBUG
            assert(f.src < output_size);
            #endif

            for (uint64_t i = 0; i < f.len; i++) out_it[pos + i] = out_it[f.src + i];
            pos += f.len;
        }
    }
}

template <typename int_t>
lz77_sss::int_text_t<int_t> lz77_sss::int_text(int_t* input, uint64_t input_size, uint16_t p)
    requires(is_int_input_v<int_t>)
{
    using sym_t = int_symbol_t<int_t>;
    sym_t* symbols = reinterpret_cast<sym_t*>(input);

    if constexpr (sizeof(sym_t) == 1) {
        return int_text_t<int_t>(symbols, input_size);
    } else {
        if (p == 0) p = omp_get_max_threads();
        uint64_t max_value = 0;

        #pragma omp parallel for num_threads(p) reduction(max : max_value)
        for (uint64_t i = 0; i < input_size; i++) max_value = std::max<uint64_t>(max_value, symbols[i]);

        return int_text_t<int_t>(symbols, input_size, max_value == UINT64_MAX ? max_value : max_value + 1);
    }
}

template <typename int_t, typename output_fnc_t>
void lz77_sss::factorize_approximate(int_t* input, uint64_t input_size, output_fnc_t output, parameters params)
    requires(is_int_input_v<int_t>)
{
    factorize_approximate(int_text(input, input_size, params.num_threads), output, params);
}

template <typename int_t, typename output_fnc_t>
void lz77_sss::factorize_exact(int_t* input, uint64_t input_size, output_fnc_t output, parameters params)
    requires(is_int_input_v<int_t>)
{
    factorize_exact(int_text(input, input_size, params.num_threads), output, params);
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::build_LCE()
{
    if (log) {
        std::cout << "building LCE data structure" << std::flush;
    }

    uint64_t lce_baseline_bytes = malloc_count_current();
    LCE = lce_t(T, tau);
    size_sss = LCE.get_sync_set().size();
    uint64_t sa_s_bytes = LCE.get_sa_s().size_in_bytes();
    uint64_t lce_bytes = malloc_count_current() - lce_baseline_bytes - sa_s_bytes;

    if (log) {
        record_phase_time("lce", time_diff_ns(time, now()));
        std::cout << " (" << format_size(lce_bytes) << ")";
        time = log_runtime(time);
        std::cout << "tau = " << tau;
        std::cout << ", SA_S size = " << format_size(sa_s_bytes) << std::endl;
        std::cout << "peak memory consumption = " << format_size(malloc_count_peak() - baseline_bytes + T.size_in_bytes()) << std::endl;
        std::cout << "input length / SSS size = " << n / (double)size_sss << std::endl;
    }
}