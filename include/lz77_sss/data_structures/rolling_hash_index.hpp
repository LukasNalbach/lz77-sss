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

#include <random>

#include <lz77_sss/lz77_sss.hpp>

template <uint8_t num_patt_lens, typename text_t>
class rolling_hash_index {
public:
    static constexpr uint32_t no_slot = std::numeric_limits<uint32_t>::max();
    static constexpr uint64_t min_bytes = 1 << 20;
    static constexpr uint64_t max_bytes = 1 << 30;
    static constexpr double min_bytes_per_input_char = 0.1;

protected:
    text_t T;
    uint64_t n = 0;
    std::array<uint64_t, num_patt_lens> patt_lens;
    std::array<uint64_t, num_patt_lens> base_pows;
    std::vector<uint8_t> H;
    uint64_t capacity = 0;
    uint64_t owner_mul = 0;
    uint64_t base = 0;
    uint8_t entry_bytes = 0;

    inline uint32_t slot_of(uint64_t fp) const
    {
        uint64_t mixed = (fp ^ (fp >> 31)) * 0x9E3779B97F4A7C15ull;
        mixed ^= mixed >> 29;
        return uint32_t(((mixed >> 32) * capacity) >> 32);
    }

    inline uint64_t entry(uint64_t i) const
    {
        const uint8_t* at = H.data() + i * entry_bytes;
        uint32_t low;
        std::memcpy(&low, at, 4);
        return entry_bytes == 4 ? uint64_t(low)
            : (uint64_t(low) | (uint64_t(at[4]) << 32));
    }

    inline void set_entry(uint64_t i, uint64_t value)
    {
        uint8_t* at = H.data() + i * entry_bytes;
        const uint32_t low = uint32_t(value);
        std::memcpy(at, &low, 4);
        if (entry_bytes == 5) at[4] = uint8_t(value >> 32);
    }

public:
    rolling_hash_index() = default;

    static uint8_t entry_bytes_for(uint64_t n)
    {
        return n <= std::numeric_limits<uint32_t>::max() ? 4 : 5;
    }

    static uint64_t size_in_bytes_for(uint64_t n, int64_t target_size_in_bytes)
    {
        const int64_t bytes_per_slot = entry_bytes_for(n);
        const int64_t min_slots = std::max<uint64_t>(min_bytes, n * min_bytes_per_input_char) / bytes_per_slot;
        const int64_t max_slots = max_bytes / bytes_per_slot;
        const int64_t target_slots = std::max<int64_t>(0, target_size_in_bytes) / bytes_per_slot;
        return std::min<int64_t>(max_slots, std::max<int64_t>(min_slots, target_slots)) * bytes_per_slot;
    }

    rolling_hash_index(
        const text_t& text, uint64_t size,
        std::array<uint64_t, num_patt_lens> patt_lens,
        int64_t target_size_in_bytes, uint16_t p)
        : T(text)
        , n(size)
        , patt_lens(patt_lens)
    {
        entry_bytes = entry_bytes_for(n);
        capacity = size_in_bytes_for(n, target_size_in_bytes) / entry_bytes;
        no_init_resize(H, capacity * entry_bytes);
        advise_huge_pages(H.data(), H.size());
        lce::util::parallel_memset(H.data(), 255, H.size(), p);

        base = std::mt19937_64(std::random_device()())() | 1;

        for (uint64_t i = 0; i < num_patt_lens; i++) {
            base_pows[i] = 1;
            for (uint64_t j = 0; j < patt_lens[i]; j++) base_pows[i] *= base;
        }
    }

    inline uint64_t num_slots() const { return capacity; }

    void set_num_owners(uint64_t num_owners)
    {
        owner_mul = (num_owners << 32) / capacity;
    }

    inline uint64_t owner(uint32_t slot) const { return (uint64_t(slot) * owner_mul) >> 32; }

    inline uint64_t exchange(uint32_t slot, uint64_t pos)
    {
        const uint64_t occ = entry(slot);
        set_entry(slot, pos);
        return occ;
    }

    inline void prefetch(uint32_t slot) const { __builtin_prefetch(H.data() + uint64_t(slot) * entry_bytes, 1); }

    void fill_slots(uint64_t beg, uint64_t end, uint32_t* slots,
        uint64_t stride, std::vector<uint64_t>& prefix) const
    {
        const uint64_t lim = std::min<uint64_t>(end + patt_lens[num_patt_lens - 1], n);
        no_init_resize(prefix, lim - beg + 1);
        uint64_t fp = 0;
        prefix[0] = 0;

        for (uint64_t j = beg; j < lim; j++) {
            fp = fp * base + uint64_t(T[j]);
            prefix[j - beg + 1] = fp;
        }

        const uint64_t* pre = prefix.data();

        for (uint64_t i = 0; i < num_patt_lens; i++) {
            const uint64_t len = patt_lens[i];
            const uint64_t pow = base_pows[i];
            const uint64_t valid_end = n > len ? std::min<uint64_t>(end, n - len) : beg;
            const uint64_t num_valid = valid_end > beg ? valid_end - beg : 0;
            uint32_t* out = slots + i * stride;

            for (uint64_t j = 0; j < num_valid; j++) {
                out[j] = slot_of(pre[j + len] - pre[j] * pow);
            }

            for (uint64_t j = num_valid; j < end - beg; j++) {
                out[j] = no_slot;
            }
        }
    }

    inline uint64_t size_in_bytes() const { return sizeof(*this) + H.size(); }
};
