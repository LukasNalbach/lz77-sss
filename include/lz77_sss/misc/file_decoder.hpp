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

#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include <lz77_sss/misc/block_io.hpp>

#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

class file_decoder {
public:
    static constexpr uint64_t default_ring_bits = 26;
    static constexpr uint64_t max_pending_copies = 1024 * 1024;

    file_decoder(const std::string& file_path, uint64_t size, bool in_ram, uint16_t p,
        uint64_t ring_bits = default_ring_bits)
        : path(file_path)
        , in_ram(in_ram)
        , ring_size(uint64_t { 1 } << std::max<uint64_t>(2, ring_bits))
        , window_size(ring_size / 4)
        , mask(in_ram ? std::numeric_limits<uint64_t>::max() : ring_size - 1)
        , p(std::max<uint16_t>(1, p))
    {
        ring = aligned_text_buffer(in_ram ? size : ring_size, 0);

        #ifdef _WIN32
        out.open(file_path, std::ios::binary | std::ios::trunc);
        ok = out.good();
        readers.resize(p);
        #else
        fd = ::open(file_path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
        ok = fd >= 0;
        #endif
    }

    file_decoder(const file_decoder&) = delete;
    file_decoder& operator=(const file_decoder&) = delete;

    ~file_decoder() { finish(); }

    bool good() const { return ok; }
    uint64_t position() const { return pos; }

    void literal(char c)
    {
        ring.data()[pos++ & mask] = c;
        if (pos == window_beg + window_size) flush_window();
    }

    void literals(const char* data, uint64_t len)
    {
        while (len > 0) {
            const uint64_t take = std::min<uint64_t>(len, window_beg + window_size - pos);
            std::memcpy(ring.data() + (pos & mask), data, take);
            pos += take;
            data += take;
            len -= take;
            if (pos == window_beg + window_size) flush_window();
        }
    }

    void copy(uint64_t src, uint64_t len)
    {
        while (len > 0) {
            const uint64_t take = std::min<uint64_t>(len, window_beg + window_size - pos);
            const uint64_t history_size = ring_size - window_size;
            const uint64_t far_end = in_ram ? window_beg : (window_beg > history_size ? window_beg - history_size : 0);
            const uint64_t far = src >= far_end ? 0 : std::min<uint64_t>(take, far_end - src);
            if (far > 0) far_copies.emplace_back(copy_op { .dst = pos, .src = src, .len = far });
            if (far < take) near_copies.emplace_back(copy_op { .dst = pos + far, .src = src + far, .len = take - far });
            pos += take;
            src += take;
            len -= take;

            if (pos == window_beg + window_size) {
                flush_window();
            } else if (far_copies.size() + near_copies.size() >= max_pending_copies) {
                resolve_copies();
            }
        }
    }

    void finish()
    {
        if (finished) return;
        finished = true;
        if (pos > window_beg) flush_window();
        ring = aligned_text_buffer();

        #ifdef _WIN32
        out.close();
        readers.clear();
        #else
        if (fd >= 0) ::close(fd);
        fd = -1;
        #endif
    }

private:
    struct copy_op {
        uint64_t dst;
        uint64_t src;
        uint64_t len;
    };

    static void repeat(char* dst, uint64_t dist, uint64_t done, uint64_t len)
    {
        while (done < len) {
            const uint64_t period = done - done % dist;
            const uint64_t take = std::min<uint64_t>(len - done, period);
            std::memcpy(dst + done, dst + done - period, take);
            done += take;
        }
    }

    void copy_near(const copy_op& op)
    {
        char* dst = ring.data() + (op.dst & mask);
        const uint64_t head = std::min<uint64_t>(op.len, op.dst - op.src);

        for (uint64_t done = 0; done < head;) {
            const uint64_t at = (op.src + done) & mask;
            const uint64_t take = in_ram ? head - done : std::min<uint64_t>(head - done, ring_size - at);
            std::memcpy(dst + done, ring.data() + at, take);
            done += take;
        }

        repeat(dst, op.dst - op.src, head, op.len);
    }

    void copy_far(const copy_op& op)
    {
        char* dst = ring.data() + (op.dst & mask);

        if (in_ram) {
            std::memcpy(dst, ring.data() + op.src, op.len);
            return;
        }

        #ifdef _WIN32
        std::ifstream& reader = readers[omp_get_thread_num()];
        if (!reader.is_open()) reader.open(path, std::ios::binary);
        reader.seekg(op.src, std::ios::beg);
        reader.read(dst, op.len);
        if (!reader.good()) ok = false;
        #else
        for (uint64_t done = 0; done < op.len;) {
            const ssize_t num_bytes = ::pread(fd, dst + done, op.len - done, off_t(op.src + done));

            if (num_bytes <= 0) {
                ok = false;
                break;
            }

            done += uint64_t(num_bytes);
        }
        #endif
    }

    void resolve_copies()
    {
        if (far_copies.size() > 1 && p > 1) {
            #pragma omp parallel for num_threads(p) schedule(dynamic, 64)
            for (uint64_t i = 0; i < far_copies.size(); i++) copy_far(far_copies[i]);
        } else {
            for (const copy_op& op : far_copies) copy_far(op);
        }

        for (const copy_op& op : near_copies) copy_near(op);
        far_copies.clear();
        near_copies.clear();
    }

    void flush_window()
    {
        resolve_copies();
        const char* data = ring.data() + (window_beg & mask);
        const uint64_t len = pos - window_beg;

        #ifdef _WIN32
        out.write(data, len);
        out.flush();
        ok = ok && out.good();
        #else
        for (uint64_t done = 0; done < len;) {
            const ssize_t num_bytes = ::pwrite(fd, data + done, len - done, off_t(window_beg + done));

            if (num_bytes <= 0) {
                ok = false;
                break;
            }

            done += uint64_t(num_bytes);
        }
        #endif

        window_beg = pos;
    }

    std::string path;
    bool in_ram = false;
    uint64_t ring_size = 0;
    uint64_t window_size = 0;
    uint64_t mask = 0;
    uint16_t p = 1;
    bool ok = false;
    bool finished = false;
    uint64_t pos = 0;
    uint64_t window_beg = 0;
    aligned_text_buffer ring;
    std::vector<copy_op> far_copies;
    std::vector<copy_op> near_copies;
    #ifdef _WIN32
    std::ofstream out;
    std::vector<std::ifstream> readers;
    #else
    int fd = -1;
    #endif
};
