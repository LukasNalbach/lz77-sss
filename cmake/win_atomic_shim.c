/*
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

#include <stddef.h>
#include <string.h>
#include <intrin.h>

#define LZ77_SSS_ATOMIC_NUM_LOCKS 64u

static volatile long lz77_sss_atomic_locks[LZ77_SSS_ATOMIC_NUM_LOCKS];

static volatile long* lz77_sss_atomic_lock_for(const void* ptr)
{
    size_t h = (size_t)ptr;
    h = (h >> 4) ^ (h >> 12);
    return &lz77_sss_atomic_locks[h & (LZ77_SSS_ATOMIC_NUM_LOCKS - 1u)];
}

static void lz77_sss_atomic_acquire(volatile long* lock)
{
    while (_InterlockedExchange(lock, 1) != 0) {
        _mm_pause();
    }
}

static void lz77_sss_atomic_release(volatile long* lock)
{
    _InterlockedExchange(lock, 0);
}

void lz77_sss_atomic_load(size_t size, void* src, void* dest, int memorder)
    __asm__("__atomic_load");
int lz77_sss_atomic_compare_exchange(size_t size, void* ptr, void* expected, void* desired,
                                     int success, int failure)
    __asm__("__atomic_compare_exchange");

void lz77_sss_atomic_load(size_t size, void* src, void* dest, int memorder)
{
    (void)memorder;
    volatile long* lock = lz77_sss_atomic_lock_for(src);
    lz77_sss_atomic_acquire(lock);
    memcpy(dest, src, size);
    lz77_sss_atomic_release(lock);
}

int lz77_sss_atomic_compare_exchange(size_t size, void* ptr, void* expected, void* desired,
                                     int success, int failure)
{
    (void)success;
    (void)failure;
    volatile long* lock = lz77_sss_atomic_lock_for(ptr);
    lz77_sss_atomic_acquire(lock);
    int matched;
    if (memcmp(ptr, expected, size) == 0) {
        memcpy(ptr, desired, size);
        matched = 1;
    } else {
        memcpy(expected, ptr, size);
        matched = 0;
    }
    lz77_sss_atomic_release(lock);
    return matched;
}
