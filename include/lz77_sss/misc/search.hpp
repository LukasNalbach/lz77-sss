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

#include <lz77_sss/misc/utils.hpp>

template <typename val_t, typename idx_t, typename fnc_t>
idx_t bin_search_min_geq(val_t value, idx_t left, idx_t right, fnc_t value_at)
{
    idx_t middle;

    while (left != right) {
        middle = left + (right - left) / 2;

        if (value <= value_at(middle)) {
            right = middle;
        } else {
            left = middle + 1;
        }
    }

    return left;
}

template <typename val_t, typename idx_t, typename fnc_t>
idx_t bin_search_max_lt(val_t value, idx_t left, idx_t right, fnc_t value_at)
{
    idx_t middle;

    while (left != right) {
        middle = left + (right - left) / 2 + 1;

        if (value_at(middle) < value) {
            left = middle;
        } else {
            right = middle - 1;
        }
    }

    return left;
}

template <typename val_t, typename idx_t, typename fnc_t>
idx_t bin_search_max_geq(val_t value, idx_t left, idx_t right, fnc_t value_at)
{
    idx_t middle;

    while (left != right) {
        middle = left + (right - left) / 2 + 1;

        if (value_at(middle) >= value) {
            left = middle;
        } else {
            right = middle - 1;
        }
    }

    return left;
}

template <typename val_t, typename idx_t, typename fnc_t>
idx_t bin_search_max_leq(val_t value, idx_t left, idx_t right, fnc_t value_at)
{
    idx_t middle;

    while (left != right) {
        middle = left + (right - left) / 2 + 1;

        if (value_at(middle) <= value) {
            left = middle;
        } else {
            right = middle - 1;
        }
    }

    return left;
}

template <typename val_t, typename idx_t, direction search_dir, typename fnc_t>
idx_t exp_search_max_geq(val_t value, idx_t left, idx_t right, fnc_t value_at)
{
    if (right == left) {
        return left;
    }

    idx_t cur_step_size = 1;

    if constexpr (search_dir == LEFT) {
        right -= cur_step_size;

        while (value_at(right) < value) {
            cur_step_size *= 2;

            if (right < left + cur_step_size) {
                cur_step_size = right - left;
                right = left;
                break;
            }

            right -= cur_step_size;
        }

        return bin_search_max_geq<val_t, idx_t>(value, right, right + cur_step_size - 1, value_at);
    } else {
        left += cur_step_size;

        while (value_at(left) >= value) {
            cur_step_size *= 2;

            if (right < cur_step_size || right - cur_step_size < left) {
                cur_step_size = right - left + 1;
                left = right + 1;
                break;
            }

            left += cur_step_size;
        }

        return bin_search_max_geq<val_t, idx_t>(value, left - cur_step_size, left - 1, value_at);
    }
}
