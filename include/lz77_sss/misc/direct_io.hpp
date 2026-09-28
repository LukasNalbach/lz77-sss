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
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <istream>
#include <new>
#include <ostream>
#include <streambuf>
#include <string>
#include <system_error>
#include <vector>

#include <lz77_sss/misc/utils.hpp>

inline constexpr uint64_t direct_io_alignment = 4096;
inline constexpr uint64_t bulk_io_size = 16 * 1024 * 1024;

#ifdef _WIN32

class direct_ifstream : public std::ifstream {
public:
    direct_ifstream() = default;

    explicit direct_ifstream(const std::string& path, std::ios_base::openmode mode = std::ios_base::in)
        : std::ifstream(path, mode | std::ios_base::binary) { }

    void open(const std::string& path, std::ios_base::openmode mode = std::ios_base::in)
    {
        std::ifstream::open(path, mode | std::ios_base::binary);
    }

    direct_ifstream(const std::string& path, uint64_t)
        : std::ifstream(path, std::ios_base::in | std::ios_base::binary) { }

    void open(const std::string& path, uint64_t)
    {
        std::ifstream::open(path, std::ios_base::in | std::ios_base::binary);
    }
};

class direct_ofstream : public std::ofstream {
public:
    direct_ofstream() = default;

    explicit direct_ofstream(const std::string& path, std::ios_base::openmode mode = std::ios_base::out)
        : std::ofstream(path, mode | std::ios_base::binary) { }

    void open(const std::string& path, std::ios_base::openmode mode = std::ios_base::out)
    {
        std::ofstream::open(path, mode | std::ios_base::binary);
    }
};

#else

#include <fcntl.h>
#include <unistd.h>

inline void drop_range_from_page_cache(int fd, uint64_t off, uint64_t len)
{
    if (fd < 0 || len == 0) return;
    sync_file_range(fd, off_t(off), off_t(len),
        SYNC_FILE_RANGE_WAIT_BEFORE | SYNC_FILE_RANGE_WRITE | SYNC_FILE_RANGE_WAIT_AFTER);
    posix_fadvise(fd, off, len, POSIX_FADV_DONTNEED);
}

inline void advise_uncached(int fd)
{
    if (fd < 0) return;
    posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
    posix_fadvise(fd, 0, 0, POSIX_FADV_NOREUSE);
}

class direct_filebuf_base : public std::streambuf {
public:
    static constexpr uint64_t drop_every = 64 * 1024 * 1024;
    static constexpr uint64_t default_buffer_length = 1024 * 1024;

    bool is_open() const { return fd >= 0; }

protected:
    void attach_fd(int descriptor, uint64_t length = default_buffer_length)
    {
        fd = descriptor;
        dropped = 0;

        if (buffer != nullptr && length != buffer_length) {
            free(buffer);
            buffer = nullptr;
        }

        buffer_length = length;

        if (buffer == nullptr) {
            void* memory = nullptr;
            if (posix_memalign(&memory, direct_io_alignment, buffer_length) != 0) memory = nullptr;
            buffer = (char*) memory;
        }
    }

    ~direct_filebuf_base() override { if (buffer != nullptr) free(buffer); }

    int open_direct(const std::string& path, int flags, mode_t mode = 0644)
    {
        int descriptor = open(path.c_str(), flags | O_DIRECT, mode);
        direct = descriptor >= 0;
        if (!direct) descriptor = open(path.c_str(), flags, mode);
        return descriptor;
    }

    int fd = -1;
    uint64_t dropped = 0;
    uint64_t buffer_length = default_buffer_length;
    bool direct = false;
    char* buffer = nullptr;
};

class direct_ofilebuf : public direct_filebuf_base {
public:
    direct_ofilebuf() = default;

    explicit direct_ofilebuf(const std::string& path, bool append = false) { open_file(path, append); }

    void open_file(const std::string& path, bool append = false)
    {
        close_file();
        attach_fd(append ? open(path.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0644)
                       : open_direct(path, O_WRONLY | O_CREAT | O_TRUNC));
        if (append) direct = false;
        if (!direct) advise_uncached(fd);
        written = append && fd >= 0 ? uint64_t(std::max<off_t>(0, lseek(fd, 0, SEEK_END))) : 0;
        on_disk = written;
        dropped = written;
        in_flight_off = 0;
        in_flight_len = 0;
        setp(buffer, buffer + buffer_length);
    }

    ~direct_ofilebuf() override { close_file(); }

    direct_ofilebuf(const direct_ofilebuf&) = delete;
    direct_ofilebuf& operator=(const direct_ofilebuf&) = delete;

    void close_file()
    {
        if (fd < 0) return;
        flush_buffer();

        if (direct) {
            uint64_t rest = uint64_t(pptr() - pbase());
            written = on_disk + rest;

            if (rest > 0) {
                std::memset(buffer + rest, 0, size_t(direct_io_alignment - rest));
                write_at(buffer, direct_io_alignment, on_disk);
            }

            if (ftruncate(fd, off_t(written)) != 0) { }
        } else {
            drop_range_from_page_cache(fd, in_flight_off, in_flight_len);
            drop_range_from_page_cache(fd, dropped, written - dropped);
        }

        setp(buffer, buffer + buffer_length);
        ::close(fd);
        fd = -1;
    }

protected:
    int_type overflow(int_type c) override
    {
        if (!flush_buffer()) return traits_type::eof();
        if (c != traits_type::eof()) { *pptr() = traits_type::to_char_type(c); pbump(1); }
        return c;
    }

    int sync() override { return flush_buffer() ? 0 : -1; }

    std::streamsize xsputn(const char* s, std::streamsize n) override
    {
        if (direct || n < std::streamsize(buffer_length)) return std::streambuf::xsputn(s, n);
        if (!flush_buffer()) return 0;
        return write_all(s, uint64_t(n)) ? n : 0;
    }

private:
    bool flush_buffer()
    {
        uint64_t num_bytes = uint64_t(pptr() - pbase());
        if (num_bytes == 0) return true;

        if (!direct) {
            setp(buffer, buffer + buffer_length);
            return write_all(buffer, num_bytes);
        }

        uint64_t whole = (num_bytes / direct_io_alignment) * direct_io_alignment;
        uint64_t rest = num_bytes - whole;

        if (whole > 0 && !write_at(buffer, whole, on_disk)) return false;
        on_disk += whole;
        if (rest > 0) std::memmove(buffer, buffer + whole, size_t(rest));

        setp(buffer, buffer + buffer_length);
        pbump(int(rest));
        return true;
    }

    bool write_at(const char* data, uint64_t size, uint64_t off)
    {
        uint64_t done = 0;

        while (done < size) {
            ssize_t num_bytes = ::pwrite(fd, data + done, size_t(size - done), off_t(off + done));
            if (num_bytes <= 0) return false;
            done += uint64_t(num_bytes);
        }

        return true;
    }

    bool write_all(const char* data, uint64_t size)
    {
        uint64_t at = 0;

        while (at < size) {
            ssize_t num_bytes = ::write(fd, data + at, size - at);
            if (num_bytes <= 0) return false;
            at += uint64_t(num_bytes);
            written += uint64_t(num_bytes);
            on_disk = written;
        }

        if (written - dropped >= drop_every) {
            sync_file_range(fd, off_t(dropped), off_t(written - dropped), SYNC_FILE_RANGE_WRITE);
            drop_range_from_page_cache(fd, in_flight_off, in_flight_len);
            in_flight_off = dropped;
            in_flight_len = written - dropped;
            dropped = written;
        }

        return true;
    }

    uint64_t written = 0;
    uint64_t on_disk = 0;
    uint64_t in_flight_off = 0;
    uint64_t in_flight_len = 0;
};

class direct_ofstream : public std::ostream {
public:
    direct_ofstream() : std::ostream(&buf) { setstate(std::ios_base::failbit); }

    explicit direct_ofstream(const std::string& path,
        std::ios_base::openmode mode = std::ios_base::out) : std::ostream(&buf)
    {
        open(path, mode);
    }

    void open(const std::string& path, std::ios_base::openmode mode = std::ios_base::out)
    {
        clear();
        buf.open_file(path, (mode & std::ios_base::app) != 0);
        if (!buf.is_open()) setstate(std::ios_base::failbit);
    }

    void close() { buf.close_file(); }
    bool is_open() const { return buf.is_open(); }

    direct_ofstream(const direct_ofstream&) = delete;
    direct_ofstream& operator=(const direct_ofstream&) = delete;
    direct_ofstream(direct_ofstream&&) = delete;
    direct_ofstream& operator=(direct_ofstream&&) = delete;

private:
    direct_ofilebuf buf;
};

class direct_ifilebuf : public direct_filebuf_base {
public:
    direct_ifilebuf() = default;

    explicit direct_ifilebuf(const std::string& path,
        uint64_t length = default_buffer_length) { open_file(path, length); }
    ~direct_ifilebuf() override { close_file(); }

    direct_ifilebuf(const direct_ifilebuf&) = delete;
    direct_ifilebuf& operator=(const direct_ifilebuf&) = delete;

    void open_file(const std::string& path, uint64_t length = default_buffer_length)
    {
        close_file();
        attach_fd(open_direct(path, O_RDONLY), length);
        if (!direct) advise_uncached(fd);
        fetched_end = 0;
        setg(buffer, buffer, buffer);
    }

    void close_file()
    {
        if (fd < 0) return;
        if (!direct) posix_fadvise(fd, dropped, fetched_end - dropped, POSIX_FADV_DONTNEED);
        ::close(fd);
        fd = -1;
    }

protected:
    int_type underflow() override
    {
        if (fd < 0) return traits_type::eof();

        if (!direct) {
            ssize_t num_bytes = ::read(fd, buffer, buffer_length);
            if (num_bytes <= 0) return traits_type::eof();

            fetched_end += uint64_t(num_bytes);
            setg(buffer, buffer, buffer + num_bytes);
            drop_consumed_pages();
            return traits_type::to_int_type(buffer[0]);
        }

        uint64_t lo = (fetched_end / direct_io_alignment) * direct_io_alignment;
        uint64_t skip = fetched_end - lo;
        uint64_t got = 0;

        while (got < buffer_length) {
            ssize_t num_bytes = ::pread(fd, buffer + got, size_t(buffer_length - got), off_t(lo + got));
            if (num_bytes <= 0) break;
            got += uint64_t(num_bytes);
        }

        if (got <= skip) return traits_type::eof();

        fetched_end = lo + got;
        setg(buffer, buffer + skip, buffer + got);
        return traits_type::to_int_type(buffer[skip]);
    }

    std::streamsize xsgetn(char* s, std::streamsize n) override
    {
        std::streamsize from_buffer = std::min<std::streamsize>(n, egptr() - gptr());

        if (from_buffer > 0) {
            std::memcpy(s, gptr(), size_t(from_buffer));
            gbump(int(from_buffer));
        }

        std::streamsize left = n - from_buffer;

        if (direct || left < std::streamsize(buffer_length))
            return from_buffer + std::streambuf::xsgetn(s + from_buffer, left);

        uint64_t got = 0;

        while (got < uint64_t(left)) {
            ssize_t num_bytes = ::read(fd, s + from_buffer + got, size_t(uint64_t(left) - got));
            if (num_bytes <= 0) break;
            got += uint64_t(num_bytes);
            fetched_end += uint64_t(num_bytes);
            drop_consumed_pages();
        }

        return from_buffer + std::streamsize(got);
    }

    pos_type seekpos(pos_type pos, std::ios_base::openmode which = std::ios_base::in) override
    {
        return seekoff(off_type(pos), std::ios_base::beg, which);
    }

    pos_type seekoff(off_type off, std::ios_base::seekdir dir,
        std::ios_base::openmode = std::ios_base::in) override
    {
        if (fd < 0) return pos_type(off_type(-1));

        off_t to;

        if (dir == std::ios_base::cur) {
            to = off_type(fetched_end - uint64_t(egptr() - gptr())) + off;
        } else if (dir == std::ios_base::beg) {
            to = off_t(off);
        } else {
            off_t end = lseek(fd, 0, SEEK_END);
            if (end < 0) return pos_type(off_type(-1));
            to = end + off;
        }

        if (to < 0) return pos_type(off_type(-1));
        if (lseek(fd, to, SEEK_SET) < 0) return pos_type(off_type(-1));

        fetched_end = uint64_t(to);
        dropped = std::min(dropped, fetched_end);
        setg(buffer, buffer, buffer);
        return pos_type(off_type(to));
    }

private:
    void drop_consumed_pages()
    {
        if (fetched_end - dropped < drop_every) return;
        posix_fadvise(fd, dropped, fetched_end - dropped, POSIX_FADV_DONTNEED);
        dropped = fetched_end;
    }

    uint64_t fetched_end = 0;
};

class direct_ifstream : public std::istream {
public:
    direct_ifstream() : std::istream(&buf) { setstate(std::ios_base::failbit); }

    explicit direct_ifstream(const std::string& path,
        std::ios_base::openmode = std::ios_base::in) : std::istream(&buf) { open(path); }

    direct_ifstream(const std::string& path, uint64_t buffer_length)
        : std::istream(&buf) { open(path, buffer_length); }

    void open(const std::string& path, std::ios_base::openmode = std::ios_base::in)
    {
        clear();
        buf.open_file(path);
        if (!buf.is_open()) setstate(std::ios_base::failbit);
    }

    void open(const std::string& path, uint64_t buffer_length)
    {
        clear();
        buf.open_file(path, buffer_length);
        if (!buf.is_open()) setstate(std::ios_base::failbit);
    }

    void close() { buf.close_file(); }
    bool is_open() const { return buf.is_open(); }

    direct_ifstream(const direct_ifstream&) = delete;
    direct_ifstream& operator=(const direct_ifstream&) = delete;
    direct_ifstream(direct_ifstream&&) = delete;
    direct_ifstream& operator=(direct_ifstream&&) = delete;

private:
    direct_ifilebuf buf;
};

#endif
