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

#include <lz77_sss/data_structures/range/range_point.hpp>
#include <lz77_sss/misc/utils.hpp>
#include <vector>

class semi_dynamic_square_grid {
public:
    static constexpr bool is_decomposed() { return false; }
    static constexpr bool is_static() { return false; }
    static constexpr bool is_dynamic() { return true; }

    using point_t = range_point_t;

    static constexpr uint64_t default_win_size = 16384;

protected:
    struct __attribute__((packed)) window {
        uint40_t beg;
        uint16_t len;
    };

    uint64_t num_points = 0;
    uint64_t win_size = 0;
    uint64_t grid_width = 0;
    uint64_t pos_max_excl = 0;
    std::vector<point_t> points;
    std::vector<window> windows;

    inline uint64_t grid_index(uint64_t x_w, uint64_t y_w) const { return grid_width * y_w + x_w; }
    inline uint64_t window_index(point_t point) const { return grid_index(point.x / win_size, point.y / win_size); }

public:
    semi_dynamic_square_grid() = default;

    semi_dynamic_square_grid(
        const std::vector<point_t>& input_points,
        uint64_t pos_max_excl, [[maybe_unused]] uint16_t p = 1,
        uint64_t win_size = default_win_size)
        : win_size(bench_win_size != 0 ? bench_win_size : win_size)
        , pos_max_excl(pos_max_excl)
    {
        grid_width = div_ceil(pos_max_excl, this->win_size);
        uint64_t num_win = grid_width * grid_width;
        no_init_resize(windows, num_win);

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < num_win; i++) {
            windows[i] = { .beg = 0, .len = 0 };
        }
        uint64_t p_idx = 0;

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < input_points.size(); i++) {
            #pragma omp atomic update
            windows[window_index(input_points[i])].len++;
        }

        for (window& w : windows) {
            w.beg = p_idx;
            p_idx += w.len;
            w.len = 0;
        }

        no_init_resize(points, input_points.size());

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < input_points.size(); i++) {
            points[i] = { };
        }
    }

    inline uint64_t size() const { return num_points; }

    inline uint64_t size_in_bytes() const
    {
        return sizeof(*this) +
            windows.size() * sizeof(window) +
            points.size() * sizeof(point_t);
    }

    inline void insert(point_t point)
    {
        window& w = windows[window_index(point)];
        points[w.beg + w.len] = point;
        w.len++;
        num_points++;
    }

    std::tuple<point_t, bool> point_in_range(
        uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const
    {
        uint64_t xw_1 = x1 / win_size;
        uint64_t xw_2 = x2 / win_size;

        uint64_t yw_1 = y1 / win_size;
        uint64_t yw_2 = y2 / win_size;

        uint64_t w_idx = grid_index(xw_1, yw_1);
        uint64_t nxt_row_offs = grid_width - (xw_2 - xw_1 + 1);

        for (uint64_t y_w = yw_1; y_w <= yw_2; y_w++) {
            for (uint64_t x_w = xw_1; x_w <= xw_2; x_w++) {
                const window& w = windows[w_idx];

                if (w.len != 0) {
                    uint64_t beg = w.beg;
                    uint64_t end = w.beg + w.len;

                    for (uint64_t i = beg; i < end; i++) {
                        const point_t& point = points[i];

                        if (x1 <= point.x && point.x <= x2 &&
                            y1 <= point.y && point.y <= y2) {
                            return { point, true };
                        }
                    }
                }

                w_idx++;
            }

            w_idx += nxt_row_offs;
        }

        return { { 0, 0 }, false };
    }

    static constexpr std::string name() { return "sdsg"; }
};

