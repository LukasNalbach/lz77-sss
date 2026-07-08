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

#include <vector>
#include <tsl/sparse_map.h>

#include <lz77_sss/data_structures/decomposed_range.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_range.hpp>
#include <lz77_sss/misc/utils.hpp>

template <typename pos_t = uint32_t>
class static_weighted_square_grid : public static_weighted_range<pos_t> {
public:
    using point_t = static_weighted_range<pos_t>::point_t;

protected:
    struct __attribute__((packed)) window {
        pos_t beg;
        uint16_t len;
    };

    pos_t num_points = 0;
    pos_t win_size = 0;
    pos_t grid_width = 0;
    pos_t pos_max = 0;
    std::vector<point_t> points;
    std::vector<window> grid;

    inline pos_t grid_index(pos_t x_w, pos_t y_w) const
    {
        return grid_width * y_w + x_w;
    }

    inline pos_t window_index(point_t p) const
    {
        return grid_index(p.x / win_size, p.y / win_size);
    }

public:
    static_weighted_square_grid() = default;

    static_weighted_square_grid(
        std::vector<point_t>& points,
        pos_t pos_max, uint16_t p = 1,
        pos_t win_size = 16384)
        : win_size(win_size)
        , pos_max(pos_max)
    {
        #if defined(BENCH_RANGE_QUERIES)
        this->win_size = cur_win_size;
        #endif

        grid_width = div_ceil(pos_max, this->win_size);
        pos_t num_win = grid_width * grid_width;
        num_points = points.size();
        grid.resize(num_win, { .beg = 0, .len = 0 });
        pos_t p_idx = 0;

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < points.size(); i++) {
            #pragma omp atomic update
            grid[window_index(points[i])].len++;
        }

        for (window& w : grid) {
            w.beg = p_idx;
            p_idx += w.len;
        }

        // this can also be done in O(n) time, but requires more space
        ips4o::parallel::sort(points.begin(), points.end(),
            [&](const point_t& p1, const point_t& p2) {
                pos_t wx_1 = window_index(p1);
                pos_t wx_2 = window_index(p2);
                return wx_1 == wx_2 ? p1.weight < p2.weight : wx_1 < wx_2;
            });

        this->points = std::move(points);
    }

    inline pos_t size() const
    {
        return num_points;
    }

    inline uint64_t size_in_bytes() const
    {
        return sizeof(this) + grid.size() * sizeof(window) + points.size() * sizeof(point_t);
    }

    std::tuple<point_t, bool> lighter_point_in_range(
        pos_t weight,
        pos_t x1, pos_t x2, pos_t y1, pos_t y2) const
    {
        pos_t xw_1 = x1 / win_size;
        pos_t xw_2 = x2 / win_size;
        pos_t yw_1 = y1 / win_size;
        pos_t yw_2 = y2 / win_size;
 
        pos_t xi_1 = xw_1 + (x1 % win_size != 0);
        pos_t yi_1 = yw_1 + (y1 % win_size != 0);
        pos_t xi_2 = xw_2 + (x2 % win_size == (win_size - 1));
        pos_t yi_2 = yw_2 + (y2 % win_size == (win_size - 1));
        bool has_contained_cells = xi_1 < xi_2 && yi_1 < yi_2;
 
        if (has_contained_cells) {
            pos_t w_idx = grid_index(xi_1, yi_1);
            pos_t nxt_row_offs = grid_width - (xi_2 - xi_1);
 
            for (pos_t y_w = yi_1; y_w < yi_2; y_w++) {
                for (pos_t x_w = xi_1; x_w < xi_2; x_w++) {
                    const window& w = grid[w_idx];

                    if (w.len != 0) {
                        const point_t& p = points[w.beg];
                        if (p.weight < weight) return { p, true };
                    }
                        
                    w_idx++;
                }
 
                w_idx += nxt_row_offs;
            }
        }

        pos_t w_idx = grid_index(xw_1, yw_1);
        pos_t nxt_row_offs = grid_width - (xw_2 - xw_1 + 1);
 
        for (pos_t y_w = yw_1; y_w <= yw_2; y_w++) {
            bool y_contained = has_contained_cells && yi_1 <= y_w && y_w < yi_2;
            pos_t x_w = xw_1;
 
            while (x_w <= xw_2) {
                if (y_contained && x_w == xi_1) {
                    w_idx += xi_2 - xi_1;
                    x_w = xi_2;
                    continue;
                }
 
                const window& w = grid[w_idx];
                pos_t e = w.beg + w.len;
 
                for (pos_t i = w.beg; i < e; i++) {
                    const point_t& p = points[i];
                    if (p.weight >= weight) break;
 
                    if (x1 <= p.x && p.x <= x2 &&
                        y1 <= p.y && p.y <= y2
                    ) return { p, true };
                }
 
                w_idx++;
                x_w++;
            }
 
            w_idx += nxt_row_offs;
        }
 
        return { { 0, 0, 0 }, false };
    }

    static constexpr std::string name()
    {
        return "static weighted square grid";
    }
};

template <typename sidx_t = uint32_t>
using decomposed_static_weighted_square_grid =
    decomposed_range<static_weighted_square_grid, sidx_t>;