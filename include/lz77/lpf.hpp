/**
 * part of LukasNalbach/lz77-sss
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 * Copyright (c) Patrick Dinklage
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

class parallel_lpf_factorizer {
private:
    static constexpr uint64_t max_size_32bit = 1ULL << 31;

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
        const std::string_view& t, const pos_t* sa, const pos_t* isa,
        pos_t b, pos_t e,
        emit_function emit_literal, emit_function emit_reference)
    {
        using spos_t = std::make_signed_t<pos_t>;
        pos_t n = t.size();

        for (pos_t i = b; i < e;) {
            pos_t cur_pos = isa[i];

            spos_t psv_pos = spos_t(cur_pos) - 1;
            while (psv_pos >= 0 && sa[psv_pos] > i) --psv_pos;
            pos_t psv_lcp = psv_pos >= 0 ? lce<pos_t>(t, i, sa[psv_pos]) : 0;

            pos_t nsv_pos = cur_pos + 1;
            while (nsv_pos < n && sa[nsv_pos] > i) ++nsv_pos;
            pos_t nsv_lcp = nsv_pos < n ? lce<pos_t>(t, i, sa[nsv_pos]) : 0;

            pos_t max_lcp = std::max(psv_lcp, nsv_lcp);
            pos_t max_pos = max_lcp == psv_lcp ? pos_t(psv_pos) : nsv_pos;
            pos_t rem = e - i;
            if (max_lcp > rem) max_lcp = rem;

            if (max_lcp >= min_ref_len) {
                assert(sa[max_pos] < i);
                emit_reference(factor(i - sa[max_pos], max_lcp));
                i += max_lcp;
            } else {
                emit_literal(factor(t[i]));
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
        pos_t n = t.size();
        auto sa = std::make_unique<pos_t[]>(n);
        auto isa = std::make_unique<pos_t[]>(n);

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

        #pragma omp parallel for num_threads(p)
        for (pos_t i = 0; i < n; i++) isa[sa[i]] = i;

        if (p <= 1) {
            factorize_range<pos_t>(t, sa.get(), isa.get(), pos_t(0), n, emit_literal, emit_reference);
            return;
        }

        #pragma omp parallel for num_threads(p)
        for (uint16_t k = 0; k < p; k++) {
            pos_t b = uint64_t(n) * k / p;
            pos_t e = uint64_t(n) * (k + 1) / p;
            std::ofstream out(tmp_file + "_" + std::to_string(k), std::ios::binary);
            auto emit = [&](factor f) { out.write((char*) &f, sizeof(factor)); };
            factorize_range<pos_t>(t, sa.get(), isa.get(), b, e, emit, emit);
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
    parallel_lpf_factorizer() : min_ref_len(2) { }

    template <std::contiguous_iterator input_t>
    requires (sizeof(std::iter_value_t<input_t>) == 1)
    void factorize(
        input_t begin, const input_t& end,
        emit_function emit_literal, emit_function emit_reference,
        uint16_t p = 0, const std::string& tmp_file = "lz77_lpf_tmp")
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
        uint16_t p = 0, const std::string& tmp_file = "lz77_lpf_tmp")
    {
        auto emit = [&](factor f) { *out++ = f; };
        factorize(begin, end, emit, emit, p, tmp_file);
    }

    uint64_t min_reference_length() const { return min_ref_len; }

    void min_reference_length(uint64_t len) { min_ref_len = len; }
};

}
