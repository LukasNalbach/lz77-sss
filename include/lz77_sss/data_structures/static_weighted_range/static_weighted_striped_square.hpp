#pragma once

#include <vector>
#include <lz77_sss/misc/utils.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_range.hpp>

template<typename pos_t = uint32_t>
class static_weighted_striped_square : static_weighted_range<pos_t> {
    public:

    using point_t = static_weighted_range<pos_t>::point_t;

    protected:
    
    pos_t n = 0;
    pos_t w = 0;
    std::vector<point_t> P;

    template <bool check_x1, bool check_x2>
    inline std::tuple<point_t, bool> lighter_point_in_segment(
        pos_t weight, pos_t x,
        pos_t x1, pos_t x2, pos_t y1, pos_t y2
    ) const {
        pos_t b = x * w;
        pos_t t = b + w;
        pos_t m;

        while (b != t) {
            m = b + (t - b) / 2;

            if (P[m].y < y1) {
                b = m + 1;
            } else if (P[m].y > y2) {
                t = m;
            } else {
                break;
            }
        }

        if (b != t) {
            for (pos_t i = m; i < t; i++) {
                if (P[i].y > y2) {
                    break;
                }

                if (P[i].weight < weight &&
                    (!check_x1 || P[i].x >= x1) &&
                    (!check_x2 || P[i].x <= x2)
                ) {
                    return {P[i], true};
                }
            }

            int64_t b_ = b;
            for (int64_t i = int64_t{m} - 1; i >= b_; i--) {
                if (P[i].y < y1) {
                    break;
                }

                if (P[i].weight < weight &&
                    (!check_x1 || P[i].x >= x1) &&
                    (!check_x2 || P[i].x <= x2)
                ) {
                    return {P[i], true};
                }
            }
        }

        return {{0, 0, 0}, false};
    }

    public:

    static_weighted_striped_square() = default;

    static_weighted_striped_square(std::vector<point_t>&& points, pos_t w) : n(points.size()), w(w) {
        P = std::move(points);
        
        ips4o::sort(P.begin(), P.end(), [&](const point_t& p1, const point_t& p2){
            pos_t x1_w = p1.x / w;
            pos_t x2_w = p2.x / w;
            return x1_w < x2_w || (x1_w == x2_w && p1.y < p2.y);
        });
    }

    inline pos_t size() const override {
        return n;
    }

    std::tuple<point_t, bool> lighter_point_in_range(pos_t weight, pos_t x1, pos_t x2, pos_t y1, pos_t y2) const override {
        pos_t x1_w = x1 / w;
        pos_t x2_w = x2 / w;

        for (pos_t x = x1_w + 1; x < x2_w; x++) {
            auto [p, result] = lighter_point_in_segment<false, false>(weight, x, x1, x2, y1, y2);

            if (result) {
                return {p, true};
            }
        }

        if (x1_w == x2_w) {
            auto [p, result] = lighter_point_in_segment<true, true>(weight, x1_w, x1, x2, y1, y2);

            if (result) {
                return {p, true};
            }
        } else {
            auto [p1, result1] = lighter_point_in_segment<true, false>(weight, x1_w, x1, x2, y1, y2);

            if (result1) {
                return {p1, true};
            }

            auto [p2, result2] = lighter_point_in_segment<false, true>(weight, x2_w, x1, x2, y1, y2);

            if (result2) {
                return {p2, true};
            }
        }

        return {{0, 0, 0}, false};
    }

    inline uint64_t size_in_bytes() const override {
        return sizeof(this) + P.size() * sizeof(point_t);
    }

    static constexpr std::string name() {
        return "static weighted striped square";
    }
};