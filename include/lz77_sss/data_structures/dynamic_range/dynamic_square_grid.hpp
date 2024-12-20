#pragma once

#include <lz77_sss/data_structures/decomposed_range.hpp>
#include <lz77_sss/data_structures/dynamic_range/dynamic_range.hpp>
#include <lz77_sss/misc/utils.hpp>
#include <vector>

template <typename pos_t = uint32_t>
class dynamic_square_grid : public dynamic_range<pos_t> {
public:
    using point_t = dynamic_range<pos_t>::point_t;

protected:
    pos_t num_points = 0;
    pos_t win_size = 0;
    pos_t grid_width = 0;
    std::vector<std::vector<point_t>> grid;

    inline pos_t grid_index(pos_t x_w, pos_t y_w) const
    {
        return grid_width * y_w + x_w;
    }

    inline pos_t window_index(point_t p) const
    {
        return grid_index(p.x / win_size, p.y / win_size);
    }

    inline std::vector<point_t>& window(point_t p)
    {
        return grid[window_index(p)];
    }

public:
    dynamic_square_grid() = default;

    dynamic_square_grid(
        std::vector<point_t>& points,
        pos_t pos_max, uint16_t p = 1,
        pos_t win_size = 16384)
        : win_size(win_size)
    {
        #if defined(BENCH_RANGE_QUERIES)
        this->win_size = cur_win_size;
        #endif

        grid_width = div_ceil(pos_max, this->win_size);
        grid.resize(grid_width * grid_width);
    }

    inline pos_t size() const override
    {
        return num_points;
    }

    inline uint64_t size_in_bytes() const override
    {
        return sizeof(this) +
            num_points * sizeof(point_t) +
            grid.size() * sizeof(std::vector<point_t>);
    }

    inline void insert(point_t p) override
    {
        window(p).emplace_back(p);
        num_points++;
    }

    std::tuple<point_t, bool> point_in_range(
        pos_t x1, pos_t x2, pos_t y1, pos_t y2) const override
    {
        pos_t xw_1 = x1 / win_size;
        pos_t xw_2 = x2 / win_size;

        pos_t yw_1 = y1 / win_size;
        pos_t yw_2 = y2 / win_size;

        pos_t w_idx = grid_index(xw_1, yw_1);
        pos_t nxt_row_offs = grid_width - (xw_2 - xw_1 + 1);

        for (pos_t y_w = yw_1; y_w <= yw_2; y_w++) {
            for (pos_t x_w = xw_1; x_w <= xw_2; x_w++) {
                for (const point_t& p : grid[w_idx]) {
                    if (x1 <= p.x && p.x <= x2 &&
                        y1 <= p.y && p.y <= y2) {
                        return { p, true };
                    }
                }

                w_idx++;
            }

            w_idx += nxt_row_offs;
        }

        return { { 0, 0 }, false };
    }

    static constexpr std::string name()
    {
        return "dynamic square grid";
    }
};

template <typename sidx_t = uint32_t>
using decomposed_dynamic_square_grid =
    decomposed_range<dynamic_square_grid, sidx_t>;