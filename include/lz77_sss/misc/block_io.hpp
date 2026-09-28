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
#include <fstream>
#include <new>
#include <string>
#include <vector>

#include <lz77_sss/misc/direct_io.hpp>
#include <lz77_sss/misc/utils.hpp>

class aligned_text_buffer {
public:
    aligned_text_buffer() = default;

    aligned_text_buffer(uint64_t size, uint64_t excess)
        : len(size)
        , capacity(((size + excess + direct_io_alignment - 1) / direct_io_alignment) * direct_io_alignment)
    {
        ptr = static_cast<char*>(::operator new[](capacity, std::align_val_t(direct_io_alignment)));
        advise_huge_pages(ptr, capacity);
        std::memset(ptr + size, 0, capacity - size);
    }

    aligned_text_buffer(aligned_text_buffer&& other) noexcept { swap(other); }

    aligned_text_buffer& operator=(aligned_text_buffer&& other) noexcept
    {
        aligned_text_buffer tmp(std::move(other));
        swap(tmp);
        return *this;
    }

    aligned_text_buffer(const aligned_text_buffer&) = delete;
    aligned_text_buffer& operator=(const aligned_text_buffer&) = delete;

    ~aligned_text_buffer()
    {
        if (ptr != nullptr) ::operator delete[](ptr, std::align_val_t(direct_io_alignment));
    }

    char* data() const { return ptr; }
    uint64_t size() const { return len; }
    bool allocated() const { return ptr != nullptr; }

private:
    void swap(aligned_text_buffer& other) noexcept
    {
        std::swap(ptr, other.ptr);
        std::swap(len, other.len);
        std::swap(capacity, other.capacity);
    }

    char* ptr = nullptr;
    uint64_t len = 0;
    uint64_t capacity = 0;
};

class positional_reader {
public:
    positional_reader(const std::string& file_path, uint16_t p)
    {
        #ifdef _WIN32
        path = file_path;
        streams.resize(std::max<uint16_t>(1, p));
        ok = std::ifstream(file_path, std::ios::binary).good();
        #else
        fd_direct = ::open(file_path.c_str(), O_RDONLY | O_DIRECT);
        fd_plain = ::open(file_path.c_str(), O_RDONLY);
        ok = fd_plain >= 0;
        #endif
    }

    positional_reader(const positional_reader&) = delete;
    positional_reader& operator=(const positional_reader&) = delete;

    ~positional_reader()
    {
        #ifndef _WIN32
        if (fd_direct >= 0) ::close(fd_direct);
        if (fd_plain >= 0) ::close(fd_plain);
        #endif
    }

    bool good() const { return ok; }

    bool read(char* buf, uint64_t len, uint64_t off, uint16_t t)
    {
        #ifdef _WIN32
        std::ifstream& in = streams[t];
        if (!in.is_open()) in.open(path, std::ios::binary);
        in.seekg(std::streamoff(off), std::ios::beg);
        in.read(buf, std::streamsize(len));
        return in.good();
        #else
        constexpr uint64_t alignment = direct_io_alignment;
        uint64_t got = 0;

        if (fd_direct >= 0 && uint64_t(uintptr_t(buf)) % alignment == 0 && off % alignment == 0) {
            const uint64_t want = ((len + alignment - 1) / alignment) * alignment;

            while (got < len) {
                const ssize_t num_bytes = ::pread(fd_direct, buf + got, size_t(want - got), off_t(off + got));

                if (num_bytes <= 0 || uint64_t(num_bytes) % alignment != 0) {
                    if (num_bytes > 0) got += uint64_t(num_bytes);
                    break;
                }

                got += uint64_t(num_bytes);
            }
        }

        while (got < len) {
            const ssize_t num_bytes = ::pread(fd_plain, buf + got, size_t(len - got), off_t(off + got));
            if (num_bytes <= 0) break;
            got += uint64_t(num_bytes);
        }

        return got >= len;
        #endif
    }

private:
    bool ok = false;
    #ifdef _WIN32
    std::string path;
    std::vector<std::ifstream> streams;
    #else
    int fd_direct = -1;
    int fd_plain = -1;
    #endif
};

class positional_writer {
public:
    positional_writer(const std::string& file_path, uint16_t p)
    {
        #ifdef _WIN32
        path = file_path;
        streams.resize(std::max<uint16_t>(1, p));
        ok = std::fstream(file_path, std::ios::in | std::ios::out | std::ios::binary).good();
        #else
        fd = ::open(file_path.c_str(), O_WRONLY);
        ok = fd >= 0;
        #endif
    }

    positional_writer(const positional_writer&) = delete;
    positional_writer& operator=(const positional_writer&) = delete;

    ~positional_writer()
    {
        #ifndef _WIN32
        if (fd >= 0) ::close(fd);
        #endif
    }

    bool good() const { return ok; }

    bool write(const char* data, uint64_t len, uint64_t off, uint16_t t)
    {
        #ifdef _WIN32
        std::fstream& out = streams[t];
        if (!out.is_open()) out.open(path, std::ios::in | std::ios::out | std::ios::binary);
        out.seekp(std::streamoff(off), std::ios::beg);
        out.write(data, std::streamsize(len));
        out.flush();
        return out.good();
        #else
        for (uint64_t done = 0; done < len;) {
            const ssize_t num_bytes = ::pwrite(fd, data + done, size_t(len - done), off_t(off + done));
            if (num_bytes <= 0) return false;
            done += uint64_t(num_bytes);
        }

        return true;
        #endif
    }

private:
    bool ok = false;
    #ifdef _WIN32
    std::string path;
    std::vector<std::fstream> streams;
    #else
    int fd = -1;
    #endif
};

template <typename fnc_t>
inline bool read_file_parallel(const std::string& path, uint64_t file_size, char* into, uint16_t p, fnc_t process)
{
    constexpr uint64_t block_size = bulk_io_size;
    const uint64_t num_blocks = (file_size + block_size - 1) / block_size;
    positional_reader reader(path, p);
    if (!reader.good()) return false;
    std::vector<aligned_text_buffer> scratch(into == nullptr ? p : 0);
    std::atomic<bool> ok = true;

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t b = 0; b < num_blocks; b++) {
        const uint16_t t = omp_get_thread_num();
        const uint64_t off = b * block_size;
        const uint64_t len = std::min<uint64_t>(block_size, file_size - off);
        if (into == nullptr && !scratch[t].allocated()) scratch[t] = aligned_text_buffer(std::min<uint64_t>(block_size, file_size), 0);
        char* buf = into == nullptr ? scratch[t].data() : into + off;

        if (!reader.read(buf, len, off, t)) {
            ok = false;
            continue;
        }

        process(buf, len, off, t);
    }

    return ok;
}
