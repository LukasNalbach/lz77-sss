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
#include <atomic>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <lz77_sss/misc/block_io.hpp>
#include <lz77_sss/misc/direct_io.hpp>
#include <lz77_sss/misc/utils.hpp>

enum fasta_mode {
    fasta_off,
    fasta_on,
    fasta_auto
};

static constexpr double max_auto_fasta_header_share = 0.1;

template <typename error_t>
inline fasta_mode parse_fasta_mode(const std::string& name, error_t error)
{
    if (name == "off") return fasta_off;
    if (name == "on") return fasta_on;
    if (name != "auto") error("error: unknown fasta mode '" + name + "'");
    return fasta_auto;
}

struct fasta_headers {
    std::string text;
    std::vector<uint64_t> seq_pos;
    std::vector<uint64_t> text_end;
    uint64_t bytes_removed = 0;
    bool active = false;

    uint64_t size() const { return seq_pos.size(); }
    uint64_t text_begin(uint64_t k) const { return k == 0 ? 0 : text_end[k - 1]; }
    std::string_view line(uint64_t k) const { return std::string_view(text.data() + text_begin(k), text_end[k] - text_begin(k)); }

    void finish(uint64_t excess)
    {
        std::string padded;
        no_init_resize_with_excess(padded, text.size(), excess);
        std::memcpy(padded.data(), text.data(), text.size());
        text = std::move(padded);
        seq_pos.shrink_to_fit();
        text_end.shrink_to_fit();
    }
};

struct fasta_state {
    bool in_header = false;
    bool at_line_start = true;
};

class fasta_filter {
public:
    fasta_filter(fasta_headers* headers, fasta_state state = fasta_state { })
        : headers(headers)
        , in_header(state.in_header)
        , at_line_start(state.at_line_start)
    { }

    static bool is_header_start(char c) { return c == '>' || c == ';'; }

    static fasta_state state_after(const char* data, uint64_t length, fasta_state state)
    {
        if (length == 0) return state;
        const auto rbeg = std::make_reverse_iterator(data + length);
        const auto rend = std::make_reverse_iterator(data);
        const auto newline = std::find(rbeg, rend, '\n');

        if (newline == rend) {
            if (state.at_line_start && is_header_start(data[0])) state.in_header = true;
            state.at_line_start = false;
            return state;
        }

        if (newline == rbeg) return fasta_state { .in_header = false, .at_line_start = true };
        return fasta_state { .in_header = is_header_start(*newline.base()), .at_line_start = false };
    }

    uint64_t header_bytes() const { return num_header_bytes; }

    uint64_t filter(const char* data, uint64_t length, char* out, uint64_t seq_off)
    {
        uint64_t i = 0;
        uint64_t pending = 0;
        uint64_t out_len = 0;

        auto keep = [&](uint64_t end) {
            if (end <= pending) return;
            if (out + out_len != data + pending) std::memmove(out + out_len, data + pending, end - pending);
            out_len += end - pending;
        };

        while (i < length) {
            if (in_header) {
                const char* newline = (const char*) std::memchr(data + i, '\n', length - i);
                const uint64_t take = newline == nullptr
                    ? length - i : uint64_t(newline - (data + i)) + 1;
                num_header_bytes += take;
                if (headers != nullptr) headers->text.append(data + i, take);
                i += take;
                pending = i;

                if (newline != nullptr) {
                    in_header = false;
                    at_line_start = true;
                    if (headers != nullptr) headers->text_end.push_back(headers->text.size());
                }

                continue;
            }

            if (at_line_start && is_header_start(data[i])) {
                keep(i);
                pending = i;
                in_header = true;
                if (headers != nullptr) headers->seq_pos.push_back(seq_off + out_len);
                continue;
            }

            const char* newline = (const char*) std::memchr(data + i, '\n', length - i);

            if (newline == nullptr) {
                i = length;
                at_line_start = false;
                break;
            }

            i = uint64_t(newline - data) + 1;
            at_line_start = true;
        }

        keep(length);
        return out_len;
    }

private:
    fasta_headers* headers;
    uint64_t num_header_bytes = 0;
    bool in_header = false;
    bool at_line_start = true;
};

struct fasta_layout {
    uint64_t block_size = 0;
    std::vector<fasta_state> start_state;
    std::vector<uint64_t> seq_off;
};

template <typename fnc_t>
inline std::optional<fasta_layout> scan_fasta(const std::string& path, uint64_t file_size, fasta_headers& headers,
    uint64_t max_header_bytes, uint16_t p, fnc_t process, uint64_t block_size = 8 * 1024 * 1024)
{
    const uint64_t num_blocks = div_ceil<uint64_t>(file_size, block_size);
    fasta_layout layout { .block_size = block_size, .start_state = std::vector<fasta_state>(num_blocks),
        .seq_off = std::vector<uint64_t>(num_blocks + 1, 0) };
    std::vector<fasta_state> end(num_blocks);
    std::vector<std::atomic<uint8_t>> ready(num_blocks);
    std::vector<fasta_headers> local(num_blocks);
    std::vector<aligned_text_buffer> buffers(p);
    std::atomic<uint64_t> header_bytes = 0;
    std::atomic<bool> exceeded = false;
    positional_reader reader(path, p);

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t b = 0; b < num_blocks; b++) {
        const uint16_t t = omp_get_thread_num();
        const uint64_t len = std::min<uint64_t>(block_size, file_size - b * block_size);
        if (!buffers[t].allocated()) buffers[t] = aligned_text_buffer(std::min<uint64_t>(block_size, file_size), 0);
        char* buf = buffers[t].data();
        const bool skip = exceeded.load(std::memory_order_relaxed);
        if (!skip) reader.read(buf, len, b * block_size, t);
        if (b > 0) spin_until([&]() { return ready[b - 1].load(std::memory_order_acquire) != 0; });
        layout.start_state[b] = b == 0 ? fasta_state { } : end[b - 1];
        end[b] = skip ? layout.start_state[b] : fasta_filter::state_after(buf, len, layout.start_state[b]);
        ready[b].store(1, std::memory_order_release);
        if (skip) continue;
        fasta_filter filter(&local[b], layout.start_state[b]);
        layout.seq_off[b + 1] = filter.filter(buf, len, buf, 0);

        if (header_bytes.fetch_add(filter.header_bytes()) + filter.header_bytes() > max_header_bytes) {
            exceeded = true;
            continue;
        }

        process(buf, layout.seq_off[b + 1], t);
    }

    if (exceeded) return std::nullopt;

    for (uint64_t b = 0; b < num_blocks; b++) {
        layout.seq_off[b + 1] += layout.seq_off[b];
        const uint64_t text_off = headers.text.size();
        headers.text.append(local[b].text);
        for (uint64_t x : local[b].seq_pos) headers.seq_pos.push_back(layout.seq_off[b] + x);
        for (uint64_t x : local[b].text_end) headers.text_end.push_back(text_off + x);
        local[b] = fasta_headers();
    }

    if (num_blocks > 0 && end.back().in_header) headers.text_end.push_back(headers.text.size());
    headers.bytes_removed = file_size - layout.seq_off.back();
    return layout;
}

template <typename fnc_t>
inline void read_fasta_sequence(const std::string& path, uint64_t file_size, const fasta_layout& layout, uint16_t p, fnc_t process)
{
    const uint64_t num_blocks = layout.start_state.size();
    std::vector<aligned_text_buffer> buffers(p);
    positional_reader reader(path, p);

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t b = 0; b < num_blocks; b++) {
        const uint16_t t = omp_get_thread_num();
        const uint64_t len = std::min<uint64_t>(layout.block_size, file_size - b * layout.block_size);
        if (!buffers[t].allocated()) buffers[t] = aligned_text_buffer(std::min<uint64_t>(layout.block_size, file_size), 0);
        char* buf = buffers[t].data();
        reader.read(buf, len, b * layout.block_size, t);
        fasta_filter filter(nullptr, layout.start_state[b]);
        process(buf, filter.filter(buf, len, buf, 0), layout.seq_off[b], t);
    }
}

inline bool looks_like_fasta(const std::string& path, uint64_t file_size)
{
    if (file_size == 0) return false;
    const uint64_t probe_size = std::min<uint64_t>(file_size, 4 * 1024 * 1024);
    std::string probe;
    no_init_resize(probe, probe_size);
    direct_ifstream in(path);
    read_fully(in, probe.data(), probe_size);
    if (!fasta_filter::is_header_start(probe[0])) return false;
    fasta_filter filter(nullptr);
    filter.filter(probe.data(), probe_size, probe.data(), 0);
    return filter.header_bytes() <= uint64_t(max_auto_fasta_header_share * probe_size);
}
