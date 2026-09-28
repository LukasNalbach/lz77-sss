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
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <text/direct_text.hpp>
#include <text/packed_text.hpp>
#include <text/split_text.hpp>

#include <lz77_sss/misc/block_io.hpp>
#include <lz77_sss/misc/direct_io.hpp>
#include <lz77_sss/misc/fasta_reader.hpp>
#include <lz77_sss/misc/log.hpp>
#include <lz77_sss/misc/utils.hpp>

using char_histogram = std::array<uint64_t, 256>;

enum text_encoding {
    auto_encoding,
    plain_encoding,
    packed_encoding,
    split_encoding
};

static constexpr double max_auto_encoded_bits = 7;
static constexpr double min_auto_split_saving = 0.2;

template <typename error_t>
inline text_encoding parse_text_encoding(const std::string& name, error_t error)
{
    if (name == "plain") return plain_encoding;
    if (name == "packed" || name == "pack") return packed_encoding;
    if (name == "split") return split_encoding;
    if (name != "auto") error("error: unknown text encoding '" + name + "'");
    return auto_encoding;
}

template <typename text_t>
class parallel_packer {
public:
    parallel_packer(text_t& text, uint16_t p)
        : text(text)
        , thread_margins(std::max<uint16_t>(1, p))
    { }

    void pack(const char* data, uint64_t len, uint64_t pos, uint16_t t)
    {
        const uint64_t end = pos + len;
        const uint64_t beg_aligned = std::min<uint64_t>(end, div_ceil<uint64_t>(pos, chunk) * chunk);
        const uint64_t end_aligned = std::max<uint64_t>(beg_aligned, (end / chunk) * chunk);
        if (beg_aligned < end_aligned) text.pack(data + (beg_aligned - pos), end_aligned - beg_aligned, beg_aligned, 1);
        for (uint64_t i = pos; i < beg_aligned; i++) thread_margins[t].emplace_back(i, data[i - pos]);
        for (uint64_t i = end_aligned; i < end; i++) thread_margins[t].emplace_back(i, data[i - pos]);
    }

    void finish()
    {
        std::vector<std::pair<uint64_t, char>> margins;

        for (std::vector<std::pair<uint64_t, char>>& local : thread_margins) {
            margins.insert(margins.end(), local.begin(), local.end());
            local = std::vector<std::pair<uint64_t, char>>();
        }

        std::sort(margins.begin(), margins.end());
        std::array<char, chunk> buf;

        for (uint64_t i = 0; i < margins.size();) {
            const uint64_t chunk_beg = (margins[i].first / chunk) * chunk;
            uint64_t k = 0;
            while (i < margins.size() && margins[i].first < chunk_beg + chunk) buf[k++] = margins[i++].second;
            text.pack(buf.data(), k, chunk_beg, 1);
        }
    }

private:
    static constexpr uint64_t chunk = text_t::chunk_symbols;

    text_t& text;
    std::vector<std::vector<std::pair<uint64_t, char>>> thread_margins;
};

inline uint64_t alphabet_size(const char_histogram& histogram)
{
    return std::count_if(histogram.begin(), histogram.end(), [](uint64_t count) { return count != 0; });
}

inline uint8_t packed_bits_per_char(const char_histogram& histogram)
{
    const uint64_t sigma = alphabet_size(histogram);
    return uint8_t(std::max<int>(1, std::bit_width(std::max<uint64_t>(sigma, 1) - 1)));
}

inline double h0_bits_per_char(const char_histogram& histogram, uint64_t size)
{
    double bits = 0;

    for (uint64_t count : histogram) {
        if (count == 0) continue;
        const double probability = count / (double)size;
        bits -= probability * std::log2(probability);
    }

    return bits;
}

inline double split_bits_per_char(const char_histogram& histogram, uint64_t size)
{
    return 8.0 * lce::text::split_text<>::size_in_bytes_for(histogram, size) / std::max<uint64_t>(1, size);
}

inline text_encoding auto_encoding_for(const char_histogram& histogram, uint64_t size)
{
    const double packed = packed_bits_per_char(histogram);
    const double split = split_bits_per_char(histogram, size);
    if (packed >= max_auto_encoded_bits && split >= max_auto_encoded_bits) return plain_encoding;
    if (packed < max_auto_encoded_bits && split > (1 - min_auto_split_saving) * packed) return packed_encoding;
    return split_encoding;
}

template <typename text_t>
inline text_t read_encoded_file(const std::string& path, uint64_t file_size, uint64_t text_size,
    const char_histogram& histogram, const fasta_layout* layout, uint16_t p = 1)
{
    text_t text(histogram, text_size);

    if (layout != nullptr) {
        parallel_packer<text_t> packer(text, p);
        read_fasta_sequence(path, file_size, *layout, p,
            [&](const char* data, uint64_t len, uint64_t pos, uint16_t t) { packer.pack(data, len, pos, t); });
        packer.finish();
    } else {
        read_file_parallel(path, text_size, nullptr, p, [&](const char* data, uint64_t len, uint64_t pos, uint16_t) {
            text.pack(data, len, pos, 1);
        });
    }

    if constexpr (requires { text.finish(p); }) text.finish(p);
    return text;
}

inline text_encoding estimate_auto_encoding(const std::string& path, uint64_t file_size, uint16_t p)
{
    const uint64_t probe_size = std::min<uint64_t>(file_size, 4 * 1024 * 1024);
    if (probe_size == 0) return auto_encoding;

    std::string probe;
    no_init_resize(probe, probe_size);
    direct_ifstream in(path);
    read_fully(in, probe.data(), probe_size);

    char_histogram histogram { };
    lce::text::packed_text::add_char_counts(histogram, probe.data(), probe_size, p);
    const double scale = double(file_size) / probe_size;
    for (uint64_t& count : histogram) count = uint64_t(std::ceil(count * scale));
    return auto_encoding_for(histogram, file_size);
}

template <typename fnc_t>
inline void with_text_from_file(const std::string& path, uint64_t file_size, text_encoding encoding,
    fasta_mode fasta, fasta_headers& headers, uint64_t excess, uint16_t p, bool log, fnc_t use,
    char_histogram* text_histogram = nullptr)
{
    auto time = now();
    char_histogram histogram { };
    aligned_text_buffer data;
    uint64_t text_size = file_size;
    std::optional<fasta_layout> layout;
    headers = fasta_headers();
    headers.active = fasta == fasta_on || (fasta == fasta_auto && looks_like_fasta(path, file_size));
    const bool scan_logged = log && headers.active;

    if (headers.active) {
        if (log) std::cout << "scanning input (" << format_size(file_size) << ")" << std::flush;
        const uint64_t max_header_bytes = fasta == fasta_on ? std::numeric_limits<uint64_t>::max()
            : uint64_t(max_auto_fasta_header_share * file_size);

        std::vector<char_histogram> partial(p, char_histogram { });

        layout = scan_fasta(path, file_size, headers, max_header_bytes, p,
            [&](const char* block, uint64_t len, uint16_t t) {
                char_histogram& local = partial[t];
                for (uint64_t i = 0; i < len; i++) local[uint8_t(block[i])]++;
            });

        if (layout) {
            text_size = layout->seq_off.back();

            for (const char_histogram& local : partial) {
                for (uint16_t c = 0; c < 256; c++) histogram[c] += local[c];
            }
        } else {
            headers = fasta_headers();
            histogram = char_histogram { };
        }
    }

    const bool scan_reads_text = !headers.active && encoding != packed_encoding && encoding != split_encoding &&
        (encoding == plain_encoding || estimate_auto_encoding(path, file_size, p) == plain_encoding);

    if (log && !headers.active && !scan_logged) {
        std::cout << (scan_reads_text ? "scanning and reading input (" : "scanning input (")
            << format_size(file_size) << ")" << std::flush;
    }

    if (scan_reads_text) {
        encoding = plain_encoding;
        data = aligned_text_buffer(text_size, excess);
        std::vector<char_histogram> partial(p, char_histogram { });

        read_file_parallel(path, text_size, data.data(), p,
            [&](const char* block, uint64_t len, uint64_t, uint16_t t) {
                char_histogram& local = partial[t];
                for (uint64_t i = 0; i < len; i++) local[uint8_t(block[i])]++;
            });

        for (const char_histogram& local : partial) {
            for (uint16_t c = 0; c < 256; c++) histogram[c] += local[c];
        }
    } else if (!headers.active) {
        std::vector<char_histogram> partial(p, char_histogram { });

        read_file_parallel(path, text_size, nullptr, p,
            [&](const char* block, uint64_t len, uint64_t, uint16_t t) {
                char_histogram& local = partial[t];
                for (uint64_t i = 0; i < len; i++) local[uint8_t(block[i])]++;
            });

        for (const char_histogram& local : partial) {
            for (uint16_t c = 0; c < 256; c++) histogram[c] += local[c];
        }
    }

    if (text_histogram != nullptr) *text_histogram = histogram;
    const uint64_t sigma = alphabet_size(histogram);
    if (headers.active) headers.finish(excess);

    if (log) {
        record_phase_time("scan_text", time_diff_ns(time, now()));
        time = log_runtime(time);
        std::cout << "sigma = " << sigma;
        std::cout << ", H0 = " << h0_bits_per_char(histogram, text_size) << " bits/char" << std::endl;

        if (headers.active) {
            std::cout << "fasta headers = " << headers.size();
            std::cout << ", header bytes = " << format_size(headers.bytes_removed) << std::endl;
        }
    }

    if (encoding == auto_encoding) {
        encoding = auto_encoding_for(histogram, text_size);
    }

    auto log_encoding = [&](double bits, uint64_t bytes) {
        if (!log) return;
        record_phase_time("encode_text", time_diff_ns(time, now()));
        std::cout << " (" << std::round(100 * bits) / 100 << " bits/char, " << std::round(1000 * bits / 8) / 10 << " %, ";
        std::cout << format_size(bytes) << ")";
        time = log_runtime(time);
    };

    if (encoding == split_encoding) {
        if (log) std::cout << "split-encoding" << std::flush;
        auto text = read_encoded_file<lce::text::split_text<>>(path, file_size, text_size, histogram,
            layout ? &*layout : nullptr, p);
        log_encoding(8.0 * text.size_in_bytes() / std::max<uint64_t>(1, text_size), text.size_in_bytes());
        use(std::move(text));
    } else if (encoding == packed_encoding) {
        if (log) std::cout << "packing" << std::flush;
        auto text = read_encoded_file<lce::text::packed_text>(path, file_size, text_size, histogram,
            layout ? &*layout : nullptr, p);
        log_encoding(text.width(), text.size_in_bytes());
        use(std::move(text));
    } else {
        if (!data.allocated()) {
            if (log) std::cout << "reading" << std::flush;
            data = aligned_text_buffer(text_size, excess);

            if (layout) {
                read_fasta_sequence(path, file_size, *layout, p, [&](const char* block, uint64_t len, uint64_t pos, uint16_t) {
                    std::memcpy(data.data() + pos, block, len);
                });
            } else {
                read_file_parallel(path, text_size, data.data(), p, [](const char*, uint64_t, uint64_t, uint16_t) { });
            }

            log_encoding(8, text_size);
        }

        use(lce::text::direct_text<char>(data.data(), text_size));
    }
}
