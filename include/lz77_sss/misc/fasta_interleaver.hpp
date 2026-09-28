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

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <lz77_sss/lz77_sss.hpp>
#include <lz77_sss/misc/fasta_reader.hpp>

template <typename out_t>
class fasta_interleaver {
public:
    template <typename factorize_t>
    fasta_interleaver(fasta_headers&& headers, out_t& out, factorize_t factorize, bool log = false)
        : out(out)
        , header_seq_pos(std::move(headers.seq_pos))
        , header_text_end(std::move(headers.text_end))
    {
        if (!header_seq_pos.empty()) {
            auto time = now();
            if (log) std::cout << "factorizing fasta headers" << std::flush;
            build_header_factors(headers.text, factorize);

            if (log) {
                std::cout << " (" << header_factors.size() << " factors)";
                log_runtime(time);
            }
        }

        headers = fasta_headers { };
    }

    void add(lz77_sss::factor f)
    {
        if (header_seq_pos.empty()) {
            emit(f);
            return;
        }

        if (f.is_literal()) {
            flush_headers();
            emit(f);
            seq_pos++;
            return;
        }

        uint64_t left = f.len;
        uint64_t src = f.src;

        while (left > 0) {
            flush_headers();
            const uint64_t headers_before_src = std::upper_bound(header_seq_pos.begin(), header_seq_pos.end(), src) - header_seq_pos.begin();
            const uint64_t take = std::min({ left, header_pos(next_header) - seq_pos, header_pos(headers_before_src) - src });
            emit(lz77_sss::factor { .src = src + (headers_before_src == 0 ? 0 : header_text_end[headers_before_src - 1]), .len = take });
            seq_pos += take;
            src += take;
            left -= take;
        }
    }

    void finish() { flush_headers(); }

    uint64_t num_headers() const { return header_seq_pos.size(); }

    uint64_t num_factors() const { return num_emitted; }

private:
    void emit(const lz77_sss::factor& f)
    {
        out(f);
        num_emitted++;
    }

    template <typename factorize_t>
    void build_header_factors(std::string& text, factorize_t factorize)
    {
        uint64_t pos = 0;
        uint64_t header = 0;

        factorize(text.data(), header_text_end.back(), [&](lz77_sss::factor f) {
            if (f.is_literal()) {
                header_factors.emplace_back(f);
                pos++;
                return;
            }

            uint64_t left = f.len;
            uint64_t src = f.src;

            while (left > 0) {
                while (header_text_end[header] <= pos) header++;
                const uint64_t src_header = std::upper_bound(header_text_end.begin(), header_text_end.end(), src) - header_text_end.begin();
                const uint64_t take = std::min({ left, header_text_end[header] - pos, header_text_end[src_header] - src });
                header_factors.emplace_back(lz77_sss::factor { .src = header_seq_pos[src_header] + src, .len = take });
                pos += take;
                src += take;
                left -= take;
            }
        });

        header_factors.shrink_to_fit();
    }

    uint64_t header_pos(uint64_t header) const
    {
        return header < header_seq_pos.size() ? header_seq_pos[header] : std::numeric_limits<uint64_t>::max();
    }

    void flush_headers()
    {
        while (next_header < header_seq_pos.size() && header_seq_pos[next_header] == seq_pos) {
            while (header_text_pos < header_text_end[next_header]) {
                const lz77_sss::factor f = header_factors[next_factor++];
                emit(f);
                header_text_pos += f.text_len();
            }

            next_header++;
        }
    }

    out_t& out;
    std::vector<uint64_t> header_seq_pos;
    std::vector<uint64_t> header_text_end;
    std::vector<lz77_sss::factor> header_factors;
    uint64_t seq_pos = 0;
    uint64_t header_text_pos = 0;
    uint64_t next_header = 0;
    uint64_t next_factor = 0;
    uint64_t num_emitted = 0;
};
