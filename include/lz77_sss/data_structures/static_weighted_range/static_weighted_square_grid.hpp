#pragma once

#include <vector>
#include <tsl/sparse_map.h>
#include <lz77_sss/misc/utils.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_range.hpp>
#include <lz77_sss/data_structures/decomposed_range.hpp>

template<typename pos_t = uint32_t>
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

    inline pos_t grid_index(pos_t x_w, pos_t y_w) const {
        return grid_width * y_w + x_w;
    }

    inline pos_t window_index(point_t p) const {
        return grid_index(p.x / win_size, p.y / win_size);
    }

    public:

    static_weighted_square_grid() = default;

    static_weighted_square_grid(
        std::vector<point_t>& points,
        pos_t pos_max, double s = 1.0
    ) : pos_max(pos_max) {
        win_size = std::ceil(std::sqrt(pos_max) / std::sqrt(s));
        grid_width = std::ceil(pos_max / (double) win_size);
        pos_t num_win = grid_width * grid_width;
        num_points = points.size();
        grid.resize(num_win, {.len = 0});
        pos_t p_idx = 0;

        for (const point_t& p : points) {
            grid[window_index(p)].len++;
        }

        for (window& w : grid) {
            w.beg = p_idx;
            p_idx += w.len;
        }

        this->points = std::move(points);

        // this can also be done in O(n) time, but requires more space
        ips4o::sort(this->points.begin(), this->points.end(),
            [&](const point_t& p1, const point_t& p2){
                return window_index(p1) < window_index(p2);
        });
    }

    inline pos_t size() const override {
        return num_points;
    }
    
    inline uint64_t size_in_bytes() const override {
        return sizeof(this) +
            grid.size() * sizeof(window) + 
            points.size() * sizeof(point_t);
    }

    std::tuple<point_t, bool> lighter_point_in_range(
        pos_t weight,
        pos_t x1, pos_t x2, pos_t y1, pos_t y2
    ) const override {
        pos_t xw_1 = x1 / win_size;
        pos_t xw_2 = x2 / win_size;

        pos_t yw_1 = y1 / win_size;
        pos_t yw_2 = y2 / win_size;

        pos_t w_idx = grid_index(xw_1, yw_1);
        pos_t nxt_row_offs = grid_width - (xw_2 - xw_1 + 1);

        for (pos_t y_w = yw_1; y_w <= yw_2; y_w++) {
            for (pos_t x_w = xw_1; x_w <= xw_2; x_w++) {
                const window& w = grid[w_idx];

                if (w.len != 0) {
                    pos_t b = w.beg;
                    pos_t e = w.beg + w.len;

                    for (pos_t i = b; i < e; i++) {
                        const point_t& p = points[i];

                        if (p.weight < weight &&
                            x1 <= p.x && p.x <= x2 &&
                            y1 <= p.y && p.y <= y2
                        ) {
                            return {p, true};
                        }
                    }
                }

                w_idx++;
            }

            w_idx += nxt_row_offs;
        }

        return {{0, 0}, false};
    }

    static constexpr std::string name() {
        return "static weighted square grid";
    }
};

template <typename sidx_t = uint32_t>
using decomposed_static_weighted_square_grid =
    decomposed_range<static_weighted_square_grid, sidx_t>;