#pragma once

#include <vector>
#include <lz77_sss/misc/utils.hpp>
#include <lz77_sss/data_structures/dynamic_range/dynamic_range.hpp>

template<typename pos_t = uint32_t>
class dynamic_square_grid : dynamic_range<pos_t> {
    public:
    
    using point_t = dynamic_range<pos_t>::point_t;

    protected:
    
    pos_t n = 0;
    pos_t w = 0;
    pos_t g = 0;
    std::vector<std::vector<point_t>> grid;

    inline pos_t grid_index(pos_t x_w, pos_t y_w) const {
        return g * y_w + x_w;
    }

    inline pos_t window_index(pos_t x, pos_t y) const {
        return grid_index(x / w, y / w);
    }

    inline std::vector<point_t>& window(pos_t x, pos_t y) {
        return grid[window_index(x, y)];
    }

    inline const std::vector<point_t>& window(pos_t x, pos_t y) const {
        return grid[window_index(x, y)];
    }

    inline std::tuple<point_t, bool> report_point(
        const std::vector<point_t>& window, 
        pos_t x1, pos_t x2, pos_t y1, pos_t y2
    ) const {
        for (const point_t& p : window) {
            if (x1 <= p.x && p.x <= x2 &&
                y1 <= p.y && p.y <= y2) {
                    return {p, true};
            }
        }

        return {{0, 0}, false};
    }

    public:

    dynamic_square_grid() = default;

    dynamic_square_grid(pos_t pos_max, double s) {
        w = std::ceil(std::sqrt(pos_max) / std::sqrt(s));
        g = std::ceil(pos_max / (double) w);
        grid.resize(g * g);
    }

    inline pos_t size() const override {
        return n;
    }

    inline uint64_t size_in_bytes() const override {
        return sizeof(this) + n * sizeof(point_t) +
            grid.size() * sizeof(std::vector<point_t>);
    }

    inline void insert(point_t p) override {
        window(p.x, p.y).emplace_back(p);
        n++;
    }

    std::tuple<point_t, bool> point_in_range(
        pos_t x1, pos_t x2, pos_t y1, pos_t y2
    ) const override {
        pos_t xw_1 = x1 / w;
        pos_t xw_2 = x2 / w;

        pos_t yw_1 = y1 / w;
        pos_t yw_2 = y2 / w;

        if (yw_1 + 1 < yw_2 && xw_2 + 1 < xw_2) {
            pos_t w = grid_index(xw_1 + 1, yw_1 + 1);
            pos_t x_r = xw_2 - xw_1 - 1;

            for (pos_t y_w = yw_1 + 1; y_w < yw_2; y_w++) {
                for (pos_t x_w = xw_1 + 1; x_w < xw_2; x_w++) {
                    for (const point_t& p : grid[w]) {
                        return {p, true};
                    }

                    w++;
                }

                w += g;
                w -= x_r;
            }
        }

        for (pos_t x = xw_1; x <= xw_2; x++) {
            auto [p, result] = report_point(
                grid[grid_index(x, yw_1)],
                x1, x2, y1, y2);

            if (result) {
                return {p, true};
            }

            if (yw_2 != yw_1) {
                auto [p2, result2] = report_point(
                    grid[grid_index(x, yw_2)],
                    x1, x2, y1, y2);

                if (result2) {
                    return {p2, true};
                }
            }            
        }

        for (pos_t y = yw_1 + 1; y < yw_2; y++) {
            auto [p, result] = report_point(
                grid[grid_index(xw_1, y)],
                x1, x2, y1, y2);

            if (result) {
                return {p, true};
            }

            if (xw_2 != xw_1) {
                auto [p2, result2] = report_point(
                    grid[grid_index(xw_2, y)],
                    x1, x2, y1, y2);

                if (result2) {
                    return {p2, true};
                }
            }            
        }

        return {{0, 0}, false};
    }

    static constexpr std::string name() {
        return "dynamic square grid";
    }
};