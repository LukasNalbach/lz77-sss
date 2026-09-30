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
#include <span>
#include <vector>
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
    static constexpr uint64_t chunks_per_thread = 8;

    template <typename sa_t, typename view_t>
    static sa_t lce(const view_t& t, sa_t i, sa_t j)
    {
        sa_t n = t.size();
        sa_t l = 0;
        while (i + l < n && j + l < n && t[i + l] == t[j + l]) ++l;
        return l;
    }

    uint64_t min_ref_len;

    template <typename sa_t, typename view_t>
    void factorize_range(
        const view_t& t, const sa_t* sa, const sa_t* isa,
        sa_t beg, sa_t end,
        emit_function emit_literal, emit_function emit_reference)
    {
        using signed_sa_t = std::make_signed_t<sa_t>;
        sa_t n = t.size();

        for (sa_t i = beg; i < end;) {
            sa_t rank = isa[i];

            signed_sa_t psv_pos = signed_sa_t(rank) - 1;
            while (psv_pos >= 0 && sa[psv_pos] > i) --psv_pos;
            sa_t psv_lcp = psv_pos >= 0 ? lce<sa_t>(t, i, sa[psv_pos]) : 0;

            sa_t nsv_pos = rank + 1;
            while (nsv_pos < n && sa[nsv_pos] > i) ++nsv_pos;
            sa_t nsv_lcp = nsv_pos < n ? lce<sa_t>(t, i, sa[nsv_pos]) : 0;

            sa_t max_lcp = std::max(psv_lcp, nsv_lcp);
            sa_t src_rank = max_lcp == psv_lcp ? sa_t(psv_pos) : nsv_pos;
            sa_t range_left = end - i;
            if (max_lcp > range_left) max_lcp = range_left;

            if (max_lcp >= min_ref_len) {
                assert(sa[src_rank] < i);
                emit_reference(factor(sa[src_rank], max_lcp));
                i += max_lcp;
            } else {
                if constexpr (sizeof(typename view_t::value_type) == 1) {
                    emit_literal(factor(t[i]));
                } else {
                    emit_literal(factor(uintmax_t(t[i]), 0));
                }

                i++;
            }
        }
    }

    template <typename sa_t, typename view_t>
    void factorize(
        const view_t& t,
        emit_function emit_literal, emit_function emit_reference,
        uint16_t p, const std::string& tmp_file_prefix)
    {
        sa_t n = t.size();
        std::vector<sa_t> sa(n);
        std::vector<sa_t> isa(n);

        if constexpr (sizeof(typename view_t::value_type) != 1) {
            using value_t = typename view_t::value_type;
            std::vector<value_t> alphabet(t.begin(), t.end());
            std::sort(alphabet.begin(), alphabet.end());
            alphabet.erase(std::unique(alphabet.begin(), alphabet.end()), alphabet.end());
            const int64_t k = std::max<int64_t>(1, alphabet.size());
            auto rank_of = [&](value_t c) { return std::lower_bound(alphabet.begin(), alphabet.end(), c) - alphabet.begin(); };

            if constexpr (std::is_same_v<sa_t, uint64_t>) {
                std::vector<int64_t> tmp(n);
                for (uint64_t i = 0; i < n; i++) tmp[i] = rank_of(t[i]);
                #ifdef LIBSAIS_OPENMP
                libsais64_long_omp(tmp.data(), (int64_t*) sa.data(), n, k, 0, p);
                #else
                libsais64_long(tmp.data(), (int64_t*) sa.data(), n, k, 0);
                #endif
            } else {
                std::vector<int32_t> tmp(n);
                for (uint64_t i = 0; i < n; i++) tmp[i] = int32_t(rank_of(t[i]));
                #ifdef LIBSAIS_OPENMP
                libsais_int_omp(tmp.data(), (int32_t*) sa.data(), n, int32_t(k), 0, p);
                #else
                libsais_int(tmp.data(), (int32_t*) sa.data(), n, int32_t(k), 0);
                #endif
            }
        } else if constexpr (std::is_same_v<sa_t, uint64_t>) {
            #ifdef LIBSAIS_OPENMP
            libsais64_omp((const uint8_t*) t.data(), (int64_t*) sa.data(), n, 0, nullptr, p);
            #else
            libsais64((const uint8_t*) t.data(), (int64_t*) sa.data(), n, 0, nullptr);
            #endif
        } else {
            #ifdef LIBSAIS_OPENMP
            libsais_omp((const uint8_t*) t.data(), (int32_t*) sa.data(), n, 0, nullptr, p);
            #else
            libsais((const uint8_t*) t.data(), (int32_t*) sa.data(), n, 0, nullptr);
            #endif
        }

        #pragma omp parallel for num_threads(p) schedule(dynamic, 65536)
        for (sa_t i = 0; i < n; i++) isa[sa[i]] = i;

        if (p <= 1) {
            factorize_range<sa_t, view_t>(t, sa.data(), isa.data(), sa_t(0), n, emit_literal, emit_reference);
            return;
        }

        const uint64_t num_ranges = std::min<uint64_t>(n, uint64_t(p) * chunks_per_thread);

        #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
        for (uint64_t k = 0; k < num_ranges; k++) {
            sa_t beg = uint64_t(n) * k / num_ranges;
            sa_t end = uint64_t(n) * (k + 1) / num_ranges;
            std::ofstream out(tmp_file_prefix + "_" + std::to_string(k), std::ios::binary);
            auto emit = [&](factor f) { out.write((char*) &f, sizeof(factor)); };
            factorize_range<sa_t, view_t>(t, sa.data(), isa.data(), beg, end, emit, emit);
        }

        for (uint64_t k = 0; k < num_ranges; k++) {
            std::string fname = tmp_file_prefix + "_" + std::to_string(k);
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
        uint16_t p = 0, const std::string& tmp_file_prefix = "lz77_lpf_tmp")
    {
        std::string_view t(begin, end);
        uint64_t n = t.size();

        #ifdef LIBSAIS_OPENMP
        if (p == 0) p = omp_get_max_threads();
        #else
        p = 1;
        #endif

        if (n < max_size_32bit) {
            factorize<uint32_t>(t, emit_literal, emit_reference, p, tmp_file_prefix);
        } else {
            factorize<uint64_t>(t, emit_literal, emit_reference, p, tmp_file_prefix);
        }
    }

    template <std::contiguous_iterator input_t>
    requires (sizeof(std::iter_value_t<input_t>) == 2 || sizeof(std::iter_value_t<input_t>) == 4)
    void factorize(
        input_t begin, const input_t& end,
        emit_function emit_literal, emit_function emit_reference,
        uint16_t p = 0, const std::string& tmp_file_prefix = "lz77_lpf_tmp")
    {
        std::span<const std::iter_value_t<input_t>> t(std::to_address(begin), uint64_t(end - begin));
        uint64_t n = t.size();

        #ifdef LIBSAIS_OPENMP
        if (p == 0) p = omp_get_max_threads();
        #else
        p = 1;
        #endif

        if (n < max_size_32bit) {
            factorize<uint32_t>(t, emit_literal, emit_reference, p, tmp_file_prefix);
        } else {
            factorize<uint64_t>(t, emit_literal, emit_reference, p, tmp_file_prefix);
        }
    }

    template <std::contiguous_iterator input_t, std::output_iterator<factor> output_t>
    requires (sizeof(std::iter_value_t<input_t>) == 1)
    void factorize(
        input_t begin, const input_t& end, output_t out,
        uint16_t p = 0, const std::string& tmp_file_prefix = "lz77_lpf_tmp")
    {
        auto emit = [&](factor f) { *out++ = f; };
        factorize(begin, end, emit, emit, p, tmp_file_prefix);
    }

    uint64_t min_reference_length() const { return min_ref_len; }

    void min_reference_length(uint64_t len) { min_ref_len = len; }
};

}
