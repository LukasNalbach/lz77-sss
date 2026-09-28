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
#include <limits>
#include <vector>

#include <lz77_sss/misc/utils.hpp>

template <typename array_t>
class min_tree {
public:
    static constexpr uint64_t block_bits = 6;
    static constexpr uint64_t block_size = uint64_t { 1 } << block_bits;

    min_tree(const array_t& data, uint64_t size, uint16_t p)
        : values(array_view(data))
        , num_values(size)
    {
        uint64_t count = size;

        while (count > 1) {
            const uint64_t blocks = (count + block_size - 1) >> block_bits;
            const std::vector<uint64_t>* below = levels.empty() ? nullptr : &levels.back();
            std::vector<uint64_t> level;
            no_init_resize(level, blocks);

            #pragma omp parallel for num_threads(p) schedule(dynamic, 256)
            for (uint64_t b = 0; b < blocks; b++) {
                const uint64_t beg = b << block_bits;
                const uint64_t end = std::min<uint64_t>(count, beg + block_size);
                uint64_t min = std::numeric_limits<uint64_t>::max();

                for (uint64_t i = beg; i < end; i++) {
                    min = std::min<uint64_t>(min, below == nullptr ? value(i) : (*below)[i]);
                }

                level[b] = min;
            }

            levels.emplace_back(std::move(level));
            count = blocks;
        }
    }

    inline uint64_t none() const { return num_values; }

    uint64_t previous_smaller(uint64_t x) const
    {
        const uint64_t v = value(x);
        const uint64_t beg = x & ~(block_size - 1);

        for (uint64_t y = x; y > beg;) {
            if (value(--y) < v) return y;
        }

        uint64_t b = x >> block_bits;

        for (uint64_t l = 0; l < levels.size(); l++) {
            const std::vector<uint64_t>& min = levels[l];
            const uint64_t group_beg = b & ~(block_size - 1);

            for (uint64_t c = b; c > group_beg;) {
                if (min[--c] < v) return last_below(l, c, v);
            }

            b >>= block_bits;
        }

        return none();
    }

    uint64_t next_smaller(uint64_t x) const
    {
        const uint64_t v = value(x);
        const uint64_t end = std::min<uint64_t>(num_values, (x | (block_size - 1)) + 1);

        for (uint64_t y = x + 1; y < end; y++) {
            if (value(y) < v) return y;
        }

        uint64_t b = x >> block_bits;

        for (uint64_t l = 0; l < levels.size(); l++) {
            const std::vector<uint64_t>& min = levels[l];
            const uint64_t group_end = std::min<uint64_t>(min.size(), (b | (block_size - 1)) + 1);

            for (uint64_t c = b + 1; c < group_end; c++) {
                if (min[c] < v) return first_below(l, c, v);
            }

            b >>= block_bits;
        }

        return none();
    }

    uint64_t size_in_bytes() const
    {
        uint64_t bytes = sizeof(*this);
        for (const std::vector<uint64_t>& level : levels) bytes += level.size() * sizeof(uint64_t);
        return bytes;
    }

private:
    inline uint64_t value(uint64_t i) const { return uint64_t(values[i]); }

    uint64_t last_below(uint64_t l, uint64_t c, uint64_t v) const
    {
        for (; l > 0; l--) {
            const std::vector<uint64_t>& min = levels[l - 1];
            uint64_t child = std::min<uint64_t>(min.size(), (c + 1) << block_bits);
            while (min[--child] >= v) { }
            c = child;
        }

        uint64_t y = std::min<uint64_t>(num_values, (c + 1) << block_bits);
        while (value(--y) >= v) { }
        return y;
    }

    uint64_t first_below(uint64_t l, uint64_t c, uint64_t v) const
    {
        for (; l > 0; l--) {
            const std::vector<uint64_t>& min = levels[l - 1];
            uint64_t child = c << block_bits;
            while (min[child] >= v) child++;
            c = child;
        }

        uint64_t y = c << block_bits;
        while (value(y) >= v) y++;
        return y;
    }

    std::vector<std::vector<uint64_t>> levels;
    array_view_t<array_t> values;
    uint64_t num_values;
};
