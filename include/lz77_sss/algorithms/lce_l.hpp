#pragma once

#include <cstdint>
#include <algorithm>
#include <lz77_sss/misc/utils.hpp>

template <typename pos_t>
pos_t lce_l_128(
    const char* T, pos_t i, pos_t j,
    pos_t lce_max = std::numeric_limits<pos_t>::max()
) {
    pos_t min_ij_p1 = std::min<pos_t>(i, j) + 1;
    lce_max = std::min<pos_t>(lce_max, min_ij_p1);

    if (i == j) [[unlikely]] {
        return lce_max;
    }

    if (i > j) [[unlikely]] {
        std::swap(i, j);
    }

    const uint128_t* ip1_ = reinterpret_cast<const uint128_t*>(T + i + 1);
    const uint128_t* min_ = ip1_ - lce_max / 16;
    const uint128_t* i_ = ip1_ - 1;
    const uint128_t* j_ = reinterpret_cast<const uint128_t*>(T + j + 1) - 1;

    while (i_ >= min_ && *i_ == *j_) {
        i_--;
        j_--;
    }

    if (i_ >= min_ || min_ij_p1 - lce_max >= 16) [[likely]] {
        return std::min<pos_t>(lce_max, 16 * ((ip1_ - i_) - 1) +
            std::countl_zero(*i_ ^ *j_) / 8);
    }

    pos_t min__ = i - (lce_max - 1);
    pos_t i__ = reinterpret_cast<const char*>(i_ + 1) - T;
    
    if (i__ <= min__) [[unlikely]] {
        return lce_max;
    }
    
    i__--;
    pos_t j__ = reinterpret_cast<const char*>(j_ + 1) - 1 - T;

    while (i__ >= min__ && T[i__] == T[j__]) {
        if (i__ == 0) [[unlikely]] {
            return lce_max;
        }

        i__--;
        j__--;
    }
    
    return i - i__;
}