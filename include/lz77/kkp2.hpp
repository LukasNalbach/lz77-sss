/**
 * part of LukasNalbach/lz77-sss
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 * Copyright (c) Patrick Dinklage
 * Copyright (c) 2013 Juha Karkkainen, Dominik Kempa and Simon J. Puglisi
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

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>

#include <libsais.h>
#include <libsais64.h>

#ifdef LIBSAIS_OPENMP
#include <omp.h>
#endif

#include <lz77/common.hpp>

namespace lz77 {

class parallel_kkp2_factorizer {
private:
    static constexpr uint64_t max_size_32bit = 1ULL << 31;

    static constexpr uint64_t stack_bits = 16;
    static constexpr uint64_t stack_size = 1ULL << stack_bits;
    static constexpr uint64_t stack_half = stack_size / 2;
    static constexpr uint64_t stack_mask = stack_size - 1;

    template <typename pos_t>
    static pos_t lce(const std::string_view& t, pos_t i, pos_t j)
    {
        pos_t n = t.length();
        pos_t l = 0;
        while (i + l < n && j + l < n && t[i + l] == t[j + l]) ++l;
        return l;
    }

    uint64_t min_ref_len;

    template <typename pos_t>
    void factorize_range(
        const std::string_view& t,
        const std::make_signed_t<pos_t>* psv_arr, const std::make_signed_t<pos_t>* nsv_arr,
        pos_t b, pos_t e,
        emit_function emit_literal, emit_function emit_reference)
    {
        using spos_t = std::make_signed_t<pos_t>;

        for (pos_t i = b + 1; i <= e;) {
            spos_t psv = psv_arr[i];
            spos_t nsv = nsv_arr[i];

            pos_t psv_lcp = psv >= 0 ? lce<pos_t>(t, i - 1, pos_t(psv) - 1) : 0;
            pos_t nsv_lcp = nsv >= 0 ? lce<pos_t>(t, i - 1, pos_t(nsv) - 1) : 0;

            pos_t max_lcp = std::max(psv_lcp, nsv_lcp);
            spos_t max_pos = max_lcp == psv_lcp ? psv : nsv;
            pos_t rem = e - (i - 1);
            if (max_lcp > rem) max_lcp = rem;

            if (max_lcp >= min_ref_len) {
                assert(max_pos >= 0);
                assert(pos_t(max_pos) < i);
                emit_reference(factor(i - max_pos, max_lcp));
                i += max_lcp;
            } else {
                emit_literal(factor(t[i - 1]));
                i++;
            }
        }
    }

    template <typename pos_t>
    void factorize(
        const std::string_view& t,
        emit_function emit_literal, emit_function emit_reference,
        uint16_t p, const std::string& tmp_file)
    {
        using spos_t = std::make_signed_t<pos_t>;

        pos_t n = t.size();
        auto cs = std::make_unique<spos_t[]>(n + 5);
        {
            auto sa = std::make_unique<pos_t[]>(n);

            if constexpr (std::is_same_v<pos_t, uint64_t>) {
                #ifdef LIBSAIS_OPENMP
                libsais64_omp((const uint8_t*) t.data(), (int64_t*) sa.get(), n, 0, nullptr, p);
                #else
                libsais64((const uint8_t*) t.data(), (int64_t*) sa.get(), n, 0, nullptr);
                #endif
            } else {
                #ifdef LIBSAIS_OPENMP
                libsais_omp((const uint8_t*) t.data(), (int32_t*) sa.get(), n, 0, nullptr, p);
                #else
                libsais((const uint8_t*) t.data(), (int32_t*) sa.get(), n, 0, nullptr);
                #endif
            }

            auto stack = std::make_unique<spos_t[]>(stack_size + 5);
            spos_t top = 0;
            stack[top] = 0;

            cs[0] = -1;
            for (pos_t i = 1; i <= n; i++) {
                spos_t sai = sa[i - 1] + 1;
                while (stack[top] > sai) --top;

                if ((top & stack_mask) == 0) {
                    if (stack[top] < 0) {
                        top = -stack[top];
                        while (top > sai) top = cs[top];
                        stack[0] = -cs[top];
                        stack[1] = top;
                        top = 1;
                    } else if (top == stack_size) {
                        for (pos_t j = stack_half; j <= stack_size; j++) {
                            stack[j - stack_half] = stack[j];
                        }
                        stack[0] = -stack[0];
                        top = stack_half;
                    }
                }

                cs[sai] = std::max(spos_t(0), stack[top]);
                ++top;
                stack[top] = sai;
            }
        }

        auto psv_arr = std::make_unique<spos_t[]>(n + 1);
        auto nsv_arr = std::make_unique<spos_t[]>(n + 1);

        cs[0] = 0;
        for (pos_t i = 1; i <= n; i++) {
            spos_t psv = cs[i];
            spos_t nsv = cs[psv];
            psv_arr[i] = psv;
            nsv_arr[i] = nsv;
            cs[i] = nsv;
            cs[psv] = i;
        }
        cs.reset();

        if (p <= 1) {
            factorize_range<pos_t>(t, psv_arr.get(), nsv_arr.get(), pos_t(0), n, emit_literal, emit_reference);
            return;
        }

        #pragma omp parallel for num_threads(p)
        for (uint16_t k = 0; k < p; k++) {
            pos_t b = uint64_t(n) * k / p;
            pos_t e = uint64_t(n) * (k + 1) / p;
            std::ofstream out(tmp_file + "_" + std::to_string(k), std::ios::binary);
            auto emit = [&](factor f) { out.write((char*) &f, sizeof(factor)); };
            factorize_range<pos_t>(t, psv_arr.get(), nsv_arr.get(), b, e, emit, emit);
        }

        for (uint16_t k = 0; k < p; k++) {
            std::string fname = tmp_file + "_" + std::to_string(k);
            std::ifstream in(fname, std::ios::binary);
            factor f;
            while (in.read((char*) &f, sizeof(factor))) {
                if (f.is_reference()) emit_reference(f);
                else emit_literal(f);
            }
            in.close();
            std::filesystem::remove(fname);
        }
    }

public:
    parallel_kkp2_factorizer() : min_ref_len(2) { }

    template <std::contiguous_iterator input_t>
    requires (sizeof(std::iter_value_t<input_t>) == 1)
    void factorize(
        input_t begin, const input_t& end,
        emit_function emit_literal, emit_function emit_reference,
        uint16_t p = 0, const std::string& tmp_file = "lz77_kkp2_tmp")
    {
        std::string_view t(begin, end);
        uint64_t n = t.size();

        #ifdef LIBSAIS_OPENMP
        if (p == 0) p = omp_get_max_threads();
        #else
        p = 1;
        #endif

        if (n < max_size_32bit) {
            factorize<uint32_t>(t, emit_literal, emit_reference, p, tmp_file);
        } else {
            factorize<uint64_t>(t, emit_literal, emit_reference, p, tmp_file);
        }
    }

    template <std::contiguous_iterator input_t, std::output_iterator<factor> output_t>
    requires (sizeof(std::iter_value_t<input_t>) == 1)
    void factorize(
        input_t begin, const input_t& end, output_t out,
        uint16_t p = 0, const std::string& tmp_file = "lz77_kkp2_tmp")
    {
        auto emit = [&](factor f) { *out++ = f; };
        factorize(begin, end, emit, emit, p, tmp_file);
    }

    uint64_t min_reference_length() const { return min_ref_len; }

    void min_reference_length(uint64_t len) { min_ref_len = len; }
};

}
