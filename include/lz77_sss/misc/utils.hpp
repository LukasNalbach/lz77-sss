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
#include <bit>
#include <cassert>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <span>
#include <string>
#include <thread>
#include <vector>
#include <random>

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

#ifndef _WIN32
#include <unistd.h>
#endif

#include <malloc_count/malloc_count.h>
#include <omp.h>
#include <rolling_hash/rolling_hash.hpp>
#include <util/bit_aligned_vector.hpp>
#include <util/memory.hpp>

using uint128_t = lce::rolling_hash::uint128_t;

using lce::util::uint40_t;
using lce::util::largest_value;
using lce::util::no_init_resize;
using lce::util::release_free_memory;
using lce::util::advise_huge_pages;
using lce::util::bit_aligned_vector;

template <typename array_t>
inline auto array_view(const array_t& arr)
{
    if constexpr (requires { arr.view(); }) {
        return arr.view();
    } else {
        return std::span(arr);
    }
}

template <typename array_t>
using array_view_t = decltype(array_view(std::declval<const array_t&>()));

inline std::chrono::steady_clock::time_point now() { return std::chrono::steady_clock::now(); }

inline void cpu_relax()
{
    #if defined(__x86_64__) || defined(_M_X64)
    _mm_pause();
    #else
    std::this_thread::yield();
    #endif
}

template <typename pred_t>
inline void spin_until(pred_t done)
{
    for (uint64_t i = 0; !done(); i++) {
        if (i < 1024) {
            cpu_relax();
        } else {
            std::this_thread::yield();
        }
    }
}

class spin_barrier {
    alignas(64) std::atomic<uint64_t> count = 0;
    alignas(64) std::atomic<uint64_t> generation = 0;
    uint64_t num = 0;

public:
    explicit spin_barrier(uint64_t num) : num(num) { }

    template <typename fnc_t>
    void wait(fnc_t on_last)
    {
        const uint64_t gen = generation.load(std::memory_order_acquire);

        if (count.fetch_add(1, std::memory_order_acq_rel) + 1 == num) {
            on_last();
            count.store(0, std::memory_order_relaxed);
            generation.store(gen + 1, std::memory_order_release);
        } else {
            spin_until([&]() { return generation.load(std::memory_order_acquire) != gen; });
        }
    }

    void wait() { wait([] { }); }
};

inline uint64_t time_diff_sec(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2)
{
    return std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
}

inline uint64_t time_diff_ns(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2)
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
}

inline uint64_t time_diff_ns(std::chrono::steady_clock::time_point t)
{
    return time_diff_ns(t, std::chrono::steady_clock::now());
}

inline void read_fully(std::istream& in, char* into, uint64_t size)
{
    uint64_t size_left = size;
    uint64_t bytes_to_read;

    while (size_left > 0) {
        bytes_to_read = std::min(size_left, (uint64_t)INT_MAX);
        in.read(into + (size - size_left), bytes_to_read);
        size_left -= bytes_to_read;
    }
}

template <typename char_t>
inline char_t uchar_to_char(uint8_t c)
{
    static_assert(sizeof(char_t) == 1);
    return *reinterpret_cast<char_t*>(&c);
}

template <typename char_t>
inline uint8_t char_to_uchar(char_t c)
{
    static_assert(sizeof(char_t) == 1);
    return *reinterpret_cast<uint8_t*>(&c);
}

template <typename uint_t, typename char_t>
inline char_t unsigned_to_char(uint_t c)
{
    static_assert(std::numeric_limits<uint_t>::min() == 0);
    static_assert(sizeof(char_t) == 1);
    return *reinterpret_cast<char_t*>(&c);
}

template <typename T, typename Alloc = std::allocator<T>>
class default_init_allocator : public Alloc {
    using Alloc::Alloc;

public:
    template <typename U>
    struct rebind { };
    template <typename U>
    void construct(U* ptr) noexcept(std::is_nothrow_default_constructible<U>::value) { ::new ((void*)(ptr)) U; }
    template <typename U, typename... Args>
    void construct(U* ptr, Args&&... args) { }
};

inline void no_init_resize(std::string& str, size_t size)
{
    (*reinterpret_cast<std::basic_string<char, std::char_traits<char>, default_init_allocator<char>>*>(&str)).resize(size);
    std::atomic_signal_fence(std::memory_order_seq_cst);
}

inline void no_init_resize_with_excess(std::string& str, size_t size, size_t excess)
{
    str.reserve(size + excess);
    no_init_resize(str, size);
    str.resize(size + excess);
    std::fill(str.end() - excess, str.end(), 0);
    str.resize(size);
}

inline std::string random_alphanumeric_string(uint64_t length)
{
    static std::string possible_chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> char_idx_distrib(0, possible_chars.size() - 1);
    
    std::string str_rand;
    str_rand.reserve(length);

    for (uint64_t i = 0; i < length; i++) {
        str_rand.push_back(possible_chars[char_idx_distrib(gen)]);
    }

    return str_rand;
}

enum direction { LEFT, RIGHT };

template <typename uint_t>
inline static uint_t div_ceil(const uint_t x, const uint_t y)
{
    if (x == 0) [[unlikely]] return 0;
    return 1 + (x - 1) / y;
}

template <typename fnc_t>
inline void parallel_chunks(uint64_t count, uint16_t p, fnc_t fnc, uint64_t chunks_per_thread = 16)
{
    if (count == 0) return;
    const uint64_t chunk = div_ceil<uint64_t>(div_ceil<uint64_t>(count, uint64_t(p) * chunks_per_thread), 64) * 64;
    const uint64_t num_chunks = div_ceil<uint64_t>(count, chunk);

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t c = 0; c < num_chunks; c++) {
        fnc(c * chunk, std::min<uint64_t>(count, (c + 1) * chunk));
    }
}

template <typename gen_t>
uint64_t random_log_uniform_size(uint64_t min_size, uint64_t max_size, gen_t& gen)
{
    std::uniform_real_distribution<double> log_distrib(
        std::log((double) std::max<uint64_t>(1, min_size)),
        std::log((double) std::max<uint64_t>(1, max_size)));
    return std::clamp<uint64_t>((uint64_t) std::llround(std::exp(log_distrib(gen))), min_size, max_size);
}

template <typename T>
static bool contains(const std::vector<T>& vec, const T& val)
{
    return std::find(vec.begin(), vec.end(), val) != vec.end();
}
