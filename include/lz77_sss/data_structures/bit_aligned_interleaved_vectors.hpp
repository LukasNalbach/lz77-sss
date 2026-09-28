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

#include <array>
#include <bit>
#include <cstdint>
#include <atomic>
#include <cstring>
#include <vector>

#include <lz77_sss/misc/utils.hpp>

enum class init_mode {
    zero,
    none
};

template <uint8_t num_fields>
class bit_aligned_interleaved_vectors {
public:
    bit_aligned_interleaved_vectors() = default;

    bit_aligned_interleaved_vectors(uint64_t size, std::array<uint64_t, num_fields> max_values, init_mode init = init_mode::zero)
        : num_entries(size)
    {
        bits_per_entry = 0;

        for (uint8_t c = 0; c < num_fields; c++) {
            field_offset[c] = bits_per_entry;
            field_width[c] = max_values[c] == 0 ? 0 : std::max<uint8_t>(1, std::bit_width(max_values[c]));
            field_mask[c] = field_width[c] == 0 ? 0 : ((uint64_t(1) << field_width[c]) - 1);
            bits_per_entry += field_width[c];
        }

        no_init_resize(words, (size * bits_per_entry) / 64 + 4);
        if (init == init_mode::zero) std::memset(words.data(), 0, words.size() * 8);
        advise_huge_pages(words.data(), words.size() * 8);
    }

    uint64_t num_words() const { return words.size(); }

    void zero_words(uint64_t from, uint64_t to) { std::memset(words.data() + from, 0, (to - from) * 8); }

    inline void init_parallel(uint64_t i, std::array<uint64_t, num_fields> values)
    {
        const uint64_t bit = i * bits_per_entry;

        if (bits_per_entry <= 64) {
            uint64_t entry = 0;

            for (uint8_t c = 0; c < num_fields; c++) {
                if (field_width[c] != 0) entry |= (values[c] & field_mask[c]) << field_offset[c];
            }

            or_bits(bit, entry);
        } else {
            for (uint8_t c = 0; c < num_fields; c++) {
                if (field_width[c] != 0) or_bits(bit + field_offset[c], values[c] & field_mask[c]);
            }
        }
    }

    template <uint8_t field_idx>
    inline uint64_t get(uint64_t i) const
    {
        if constexpr (field_idx >= num_fields) return 0;
        if (field_width[field_idx] == 0) [[unlikely]] return 0;
        const uint64_t bit = i * bits_per_entry + field_offset[field_idx];
        uint64_t block;
        std::memcpy(&block, reinterpret_cast<const uint8_t*>(words.data()) + (bit >> 3), 8);
        return (block >> (bit & 7)) & field_mask[field_idx];
    }

    inline uint64_t entry(uint64_t i) const
    {
        const uint64_t bit = i * bits_per_entry;
        const uint8_t* at = reinterpret_cast<const uint8_t*>(words.data()) + (bit >> 3);
        uint64_t lo;
        uint64_t hi;
        std::memcpy(&lo, at, 8);
        std::memcpy(&hi, at + 8, 8);
        const uint8_t shift = uint8_t(bit & 7);
        return (lo >> shift) | (hi << (63 - shift) << 1);
    }

    template <uint8_t field_idx>
    inline uint64_t field(uint64_t entry) const
    {
        if constexpr (field_idx >= num_fields) return 0;
        if (field_width[field_idx] == 0) [[unlikely]] return 0;
        return (entry >> field_offset[field_idx]) & field_mask[field_idx];
    }

    template <uint8_t field_idx>
    inline void set(uint64_t i, uint64_t value)
    {
        if constexpr (field_idx >= num_fields) return;
        if (field_width[field_idx] == 0) [[unlikely]] return;
        const uint64_t bit = i * bits_per_entry + field_offset[field_idx];
        uint8_t* at = reinterpret_cast<uint8_t*>(words.data()) + (bit >> 3);
        uint64_t block;
        std::memcpy(&block, at, 8);
        block = (block & ~(field_mask[field_idx] << (bit & 7))) | ((value & field_mask[field_idx]) << (bit & 7));
        std::memcpy(at, &block, 8);
    }

    template <uint8_t field_idx>
    inline void set_parallel(uint64_t i, uint64_t value)
    {
        if constexpr (field_idx >= num_fields) return;
        if (field_width[field_idx] == 0) [[unlikely]] return;
        const uint64_t bit = i * bits_per_entry + field_offset[field_idx];
        const uint64_t word = bit >> 6;
        const uint64_t offset = bit & 63;
        const uint64_t val = value & field_mask[field_idx];
        const uint64_t bits_lo = std::min<uint64_t>(64 - offset, field_width[field_idx]);
        const uint64_t mask_lo = (bits_lo == 64 ? ~uint64_t(0) : ((uint64_t(1) << bits_lo) - 1)) << offset;
        store_bits(words[word], mask_lo, (val << offset) & mask_lo);

        if (bits_lo < field_width[field_idx]) {
            const uint64_t mask_hi = (uint64_t(1) << (field_width[field_idx] - bits_lo)) - 1;
            store_bits(words[word + 1], mask_hi, (val >> bits_lo) & mask_hi);
        }
    }

    void reset(uint64_t capacity, std::array<uint64_t, num_fields> max_values)
    {
        bits_per_entry = 0;

        for (uint8_t c = 0; c < num_fields; c++) {
            field_offset[c] = bits_per_entry;
            field_width[c] = max_values[c] == 0 ? 0 : std::max<uint8_t>(1, std::bit_width(max_values[c]));
            field_mask[c] = field_width[c] == 0 ? 0 : ((uint64_t(1) << field_width[c]) - 1);
            bits_per_entry += field_width[c];
        }

        num_entries = 0;
        no_init_resize(words, (capacity * bits_per_entry) / 64 + 4);
        std::memset(words.data(), 0, words.size() * 8);
    }

    inline void push_back(std::array<uint64_t, num_fields> values)
    {
        if (((num_entries + 1) * bits_per_entry) / 64 + 4 > words.size()) [[unlikely]] {
            words.resize(std::max<uint64_t>(16, words.size() * 2), 0);
        }

        for (uint8_t c = 0; c < num_fields; c++) {
            if (field_width[c] == 0) continue;
            const uint64_t bit = num_entries * bits_per_entry + field_offset[c];
            uint8_t* at = reinterpret_cast<uint8_t*>(words.data()) + (bit >> 3);
            uint64_t block;
            std::memcpy(&block, at, 8);
            block = (block & ~(field_mask[c] << (bit & 7))) | ((values[c] & field_mask[c]) << (bit & 7));
            std::memcpy(at, &block, 8);
        }

        num_entries++;
    }

    void clear()
    {
        words.clear();
        words.shrink_to_fit();
        num_entries = 0;
    }

    uint64_t size() const { return num_entries; }
    uint64_t size_in_bytes() const { return words.size() * 8; }
    uint8_t entry_width() const { return bits_per_entry; }

private:
    static void store_bits(uint64_t& word, uint64_t mask, uint64_t value)
    {
        std::atomic_ref<uint64_t> ref(word);
        uint64_t old = ref.load(std::memory_order_relaxed);
        while (!ref.compare_exchange_weak(old, (old & ~mask) | value, std::memory_order_relaxed)) { }
    }

    inline void or_bits(uint64_t bit, uint64_t value)
    {
        const uint64_t word = bit >> 6;
        const uint64_t offset = bit & 63;
        std::atomic_ref<uint64_t>(words[word]).fetch_or(value << offset, std::memory_order_relaxed);

        if (offset != 0 && (value >> (64 - offset)) != 0) {
            std::atomic_ref<uint64_t>(words[word + 1]).fetch_or(value >> (64 - offset), std::memory_order_relaxed);
        }
    }

    std::vector<uint64_t> words;
    std::array<uint64_t, num_fields> field_mask {};
    std::array<uint8_t, num_fields> field_width {};
    std::array<uint8_t, num_fields> field_offset {};
    uint8_t bits_per_entry = 0;
    uint64_t num_entries = 0;
};
