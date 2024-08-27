#pragma once

#include <vector>
#include <tsl/sparse_map.h>
#include <lz77_sss/misc/utils.hpp>
#include <sdsl/sd_vector.hpp>
#include <lz77_sss/data_structures/dynamic_range/dynamic_range.hpp>

template<typename pos_t = uint32_t>
class dynamic_square_grid_sdvec : dynamic_range<pos_t> {
    public:
    
    using point_t = dynamic_range<pos_t>::point_t;

    protected:
    
    pos_t n = 0;
    pos_t w = 0;
    pos_t g = 0;
    std::vector<point_t> points;
    std::vector<pos_t> window_start;
    sdsl::sd_vector<> active_windows;
    sdsl::sd_vector<>::rank_1_type active_windows_rank_1;
    sdsl::sd_vector<>::select_1_type active_windows_select_1;

    inline pos_t next_active_window_index(pos_t active_window_index) const {
        return active_windows_select_1.select(active_window_index + 1);
    }

    inline pos_t grid_index(pos_t x_w, pos_t y_w) const {
        return g * y_w + x_w;
    }

    inline pos_t active_window_index(pos_t x, pos_t y) const {
        return active_windows_rank_1.rank(grid_index(x / w, y / w));
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

    dynamic_square_grid_sdvec() = default;

    dynamic_square_grid_sdvec(
        const std::vector<point_t>& points,
        pos_t pos_max, double s
    ) {
        w = std::ceil(std::sqrt(pos_max) / std::sqrt(s));
        g = std::ceil(pos_max / (double) w);
        pos_t num_windows = g * g;
        tsl::sparse_map<pos_t, pos_t, std::identity> freq;

        for (const point_t& p : points) {
            pos_t w_idx = window_index(p.x, p.y);
            auto it = freq.find(w_idx);

            if (it == freq.end()) {
                freq.try_emplace(w_idx, 1);
            } else {
                freq[w_idx]++;
            }
        }

        pos_t num_active_windows = freq.size();

        std::vector<std::pair<uint64_t, pos_t>> freq_vec(
            freq.begin(), freq.end()
        );

        freq.clear();

        ips4o::sort(freq_vec.begin(), freq_vec.end(), [](auto p1, auto p2){
            return p1.second < p2.second;
        });

        window_start.reserve(num_active_windows);
        sdsl::sd_vector_builder active_windows_builder(num_windows, num_active_windows);
        pos_t cur_w_start = 0;

        for (uint64_t i = 0; i < num_active_windows; i++) {
            window_start.emplace_back(cur_w_start);
            cur_w_start += freq_vec[i].second;
            active_windows_builder.set(i);
        }

        freq_vec.clear();
        freq_vec.shrink_to_fit();
        active_windows = sdsl::sd_vector<>(active_windows_builder);
        active_windows_rank_1.set_vector(&active_windows);
        active_windows_rank_1.set_vector(&active_windows);
        points.reserve(n);
        
        for (pos_t i = 0; i < n; i++) {
            points.emplace_back(point_t{.x = pos_max, .y = pos_max});
        }
    }

    inline pos_t size() const override {
        return n;
    }
    
    inline uint64_t size_in_bytes() const override {
        return
            sizeof(this) +
            sdsl::size_in_bytes(active_windows) + 
            points.size() * sizeof(point_t) +
            points.size() * sizeof(pos_t);
    }

    inline void insert(point_t p) override {
        points[window_start[active_window_index(p.x, p.y)]++] = p;
        n++;
    }

    std::tuple<point_t, bool> point_in_range(pos_t x1, pos_t x2, pos_t y1, pos_t y2) const override {
        /*
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

                    w = next_active_window_index(w);
                }

                w += g;
                w -= x_r;
            }
        }

        for (pos_t x = xw_1; x <= xw_2; x++) {
            auto [p, result] = report_point(grid[grid_index(x, yw_1)], x1, x2, y1, y2);

            if (result) {
                return {p, true};
            }

            if (yw_2 != yw_1) {
                auto [p2, result2] = report_point(grid[grid_index(x, yw_2)], x1, x2, y1, y2);

                if (result2) {
                    return {p2, true};
                }
            }            
        }

        for (pos_t y = yw_1 + 1; y < yw_2; y++) {
            auto [p, result] = report_point(grid[grid_index(xw_1, y)], x1, x2, y1, y2);

            if (result) {
                return {p, true};
            }

            if (xw_2 != xw_1) {
                auto [p2, result2] = report_point(grid[grid_index(xw_2, y)], x1, x2, y1, y2);

                if (result2) {
                    return {p2, true};
                }
            }            
        }

        return {{0, 0}, false};
        */
    }

    static constexpr std::string name() {
        return "dynamic compact square grid";
    }
};