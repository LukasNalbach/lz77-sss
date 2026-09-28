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

#include <atomic>
#include <bit>
#include <cstdint>
#include <vector>

#include <lz77_sss/data_structures/bit_aligned_interleaved_vectors.hpp>
#include <lz77_sss/misc/utils.hpp>

class interval_hash_set {
public:
    interval_hash_set() = default;

    interval_hash_set(uint64_t num_intervals, uint64_t beg_max_excl, uint64_t max_width)
        : capacity(num_intervals + num_intervals / 3 + 8)
        , slots(capacity, { beg_max_excl + 1, max_width, tag_mask })
    {
        claimed.resize(capacity / 64 + 1, 0);
        narrow = slots.entry_width() <= 64;
    }

    inline void insert_parallel(uint64_t hash, uint64_t b, uint64_t e)
    {
        const uint64_t mix = mix_of(hash);
        uint64_t i = slot_of(mix);

        while (!claim(i)) {
            if (++i == capacity) [[unlikely]] i = 0;
        }

        slots.template set_parallel<0>(i, b + 1);
        slots.template set_parallel<1>(i, e - b);
        slots.template set_parallel<2>(i, mix & tag_mask);
    }

    void finish()
    {
        claimed.clear();
        claimed.shrink_to_fit();
    }

    template <typename fnc_t>
    inline bool find(uint64_t hash, fnc_t matches, uint64_t& b, uint64_t& e) const
    {
        const uint64_t mix = mix_of(hash);
        const uint64_t tag = mix & tag_mask;
        uint64_t i = slot_of(mix);

        while (true) {
            uint64_t beg;
            uint64_t width;
            uint64_t slot_tag;
            slot_at(i, beg, width, slot_tag);
            if (beg == 0) return false;

            if (slot_tag == tag && matches(beg - 1)) {
                b = beg - 1;
                e = b + width;
                return true;
            }

            if (++i == capacity) [[unlikely]] i = 0;
        }
    }

    uint64_t size_in_bytes() const { return sizeof(*this) + slots.size_in_bytes(); }

private:
    static constexpr uint64_t tag_mask = 0xFF;

    inline static uint64_t mix_of(uint64_t hash)
    {
        return hash * 0x9E3779B97F4A7C15ull;
    }

    inline uint64_t slot_of(uint64_t mix) const
    {
        return uint64_t((uint128_t(mix) * capacity) >> 64);
    }

    inline void slot_at(uint64_t i, uint64_t& beg, uint64_t& width, uint64_t& tag) const
    {
        if (narrow) [[likely]] {
            const uint64_t slot = slots.entry(i);
            beg = slots.template field<0>(slot);
            width = slots.template field<1>(slot);
            tag = slots.template field<2>(slot);
        } else {
            beg = slots.template get<0>(i);
            width = slots.template get<1>(i);
            tag = slots.template get<2>(i);
        }
    }

    inline bool claim(uint64_t i)
    {
        const uint64_t bit = uint64_t(1) << (i & 63);
        std::atomic_ref<uint64_t> ref(claimed[i >> 6]);
        return (ref.fetch_or(bit, std::memory_order_relaxed) & bit) == 0;
    }

    uint64_t capacity = 0;
    bool narrow = false;
    bit_aligned_interleaved_vectors<3> slots;
    std::vector<uint64_t> claimed;
};
