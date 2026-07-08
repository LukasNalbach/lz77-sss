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
#include <lz77_sss/misc/utils.hpp>

template <typename pos_t, typename char_t>
pos_t lce_l_64(
    const char_t* T, pos_t i, pos_t j,
    pos_t lce_max = std::numeric_limits<pos_t>::max())
{
    pos_t min_ij_p1 = std::min<pos_t>(i, j) + 1;
    lce_max = std::min<pos_t>(lce_max, min_ij_p1);

    if (i == j) [[unlikely]] {
        return lce_max;
    }

    if (i > j) [[unlikely]] {
        std::swap(i, j);
    }

    const uint64_t* ip1_ = reinterpret_cast<const uint64_t*>(T + i + 1);
    const uint64_t* min_ = ip1_ - lce_max / sizeof(uint64_t);
    const uint64_t* i_ = ip1_ - 1;
    const uint64_t* j_ = reinterpret_cast<const uint64_t*>(T + j + 1) - 1;

    while (i_ >= min_ && *i_ == *j_) {
        i_--;
        j_--;
    }

    if (i_ >= min_ || min_ij_p1 - lce_max >= sizeof(uint64_t)) [[likely]] {
        return std::min<pos_t>(lce_max, sizeof(uint64_t) * ((ip1_ - i_) - 1) +
            std::countl_zero(*i_ ^ *j_) / 8);
    }

    pos_t min__ = i - (lce_max - 1);
    pos_t i__ = reinterpret_cast<const char_t*>(i_ + 1) - T;

    if (i__ <= min__) [[unlikely]] {
        return lce_max;
    }

    i__--;
    pos_t j__ = reinterpret_cast<const char_t*>(j_ + 1) - 1 - T;

    while (i__ >= min__ && T[i__] == T[j__]) {
        if (i__ == 0) [[unlikely]] {
            return lce_max;
        }

        i__--;
        j__--;
    }

    return i - i__;
}