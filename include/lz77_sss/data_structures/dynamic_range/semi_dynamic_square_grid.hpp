#pragma once

#include <vector>
#include <tsl/sparse_map.h>
#include <lz77_sss/misc/utils.hpp>
#include <lz77_sss/data_structures/dynamic_range/dynamic_range.hpp>

template<typename pos_t = uint32_t>
class semi_dynamic_square_grid : dynamic_range<pos_t> {
    public:
    
    using point_t = dynamic_range<pos_t>::point_t;

    protected:
    
    pos_t n = 0;
    pos_t w = 0;
    pos_t g = 0;
    pos_t pos_max = 0;
    std::vector<point_t> points;
    std::vector<pos_t> cur_w_end;

    inline pos_t grid_index(pos_t x_w, pos_t y_w) const {
        return g * y_w + x_w;
    }

    inline pos_t window_index(pos_t x, pos_t y) const {
        return grid_index(x / w, y / w);
    }

    template <
        bool check_x1, bool check_x2,
        bool check_y1, bool check_y2
    > inline std::tuple<point_t, bool> report_point(
        pos_t w_idx,
        pos_t x1, pos_t x2,
        pos_t y1, pos_t y2
    ) const {
        if (w_idx == 0 ? cur_w_end[w_idx] == 0 :
            (cur_w_end[w_idx - 1] == cur_w_end[w_idx])
        ) {
            return {{0, 0}, false};
        }

        if constexpr (
            check_x1 || check_x2 ||
            check_y1 || check_y2
        ) {
            int64_t i = cur_w_end[w_idx] - 1;
            int64_t limit = w_idx == 0 ? 0 : cur_w_end[w_idx - 1];

            while (i >= limit && points[i].x != pos_max) {
                const point_t& p = points[i];

                if ((!check_x1 || x1 <= p.x) &&
                    (!check_x2 || x2 >= p.x) &&
                    (!check_y1 || y1 <= p.y) &&
                    (!check_y2 || y2 >= p.y)
                ) {
                    return {p, true};
                }
                
                i--;
            }
        } else {
            const point_t p = points[cur_w_end[w_idx] - 1];

            if (p.x != pos_max) {
                return {p, true};
            }
        }

        return {{0, 0}, false};
    }

    public:

    semi_dynamic_square_grid() = default;

    semi_dynamic_square_grid(
        const std::vector<point_t>& points,
        pos_t pos_max, double s
    ) : pos_max(pos_max) {
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

        std::vector<std::pair<uint64_t, pos_t>> freq_vec(
            freq.begin(), freq.end());
        freq.clear();

        ips4o::sort(freq_vec.begin(), freq_vec.end(),
            [](auto p1, auto p2){return p1.second < p2.second;});

        cur_w_end.reserve(num_windows);
        pos_t last_w_idx = 0;
        pos_t cur_w_start = 0;

        for (auto& p : freq_vec) {
            while (last_w_idx <= p.first) {
                cur_w_end.emplace_back(cur_w_start);
                last_w_idx++;
            }
            
            cur_w_start += p.second;
        }

        while (last_w_idx < num_windows) {
            cur_w_end.emplace_back(points.size());
            last_w_idx++;
        }
        
        for (point_t& p : points) {
            this->points.emplace_back(
                point_t{.x = pos_max, .y = pos_max});
        }
    }

    inline pos_t size() const override {
        return n;
    }
    
    inline uint64_t size_in_bytes() const override {
        return
            sizeof(this) +
            cur_w_end.size() * sizeof(pos_t) + 
            points.size() * sizeof(point_t);
    }

    inline void insert(point_t p) override {
        points[cur_w_end[window_index(p.x, p.y)]++] = p;
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
                    auto [p, result] = report_point<
                        false, false, false, false
                        >(w, x1, x2, y1, y2);
                    if (result) return {p, true};
                    w++;
                }

                w += g;
                w -= x_r;
            }
        }

        if (xw_1 + 1 < xw_2) {
            if (yw_1 == yw_2) {
                for (pos_t x = xw_1 + 1; x < xw_2; x++) {
                    auto [p, result] = report_point<
                        false, false, true, true
                        >(grid_index(x, yw_1), x1, x2, y1, y2);
                    if (result) return {p, true};
                }
            } else {
                for (pos_t x = xw_1 + 1; x < xw_2; x++) {
                    auto [p, result] = report_point<
                        false, false, true, false
                        >(grid_index(x, yw_1), x1, x2, y1, y2);
                    if (result) return {p, true};
                }

                for (pos_t x = xw_1 + 1; x < xw_2; x++) {
                    auto [p, result] = report_point<
                        false, false, false, true
                        >(grid_index(x, yw_2), x1, x2, y1, y2);
                    if (result) return {p, true};
                }
            }
        }

        if (yw_1 + 1 < yw_2) {
            if (xw_1 == xw_2) {
                for (pos_t y = yw_1 + 1; y < yw_2; y++) {
                    auto [p, result] = report_point<
                        true, true, false, false
                        >(grid_index(xw_1, y), x1, x2, y1, y2);
                    if (result) return {p, true};
                }
            } else {
                for (pos_t y = yw_1 + 1; y < yw_2; y++) {
                    auto [p, result] = report_point<
                        true, false, false, false
                        >(grid_index(xw_1, y), x1, x2, y1, y2);
                    if (result) return {p, true};
                }

                for (pos_t y = yw_1 + 1; y < yw_2; y++) {
                    auto [p, result] = report_point<
                        false, true, false, false
                        >(grid_index(xw_2, y), x1, x2, y1, y2);
                    if (result) return {p, true};
                }
            }
        }

        return {{0, 0}, false};
    }

    static constexpr std::string name() {
        return "semi-dynamic square grid";
    }
};