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

#include <lz77_sss/data_structures/decomposed_range.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_range.hpp>
#include <lz77_sss/misc/utils.hpp>
#include <vector>

template <typename pos_t = uint32_t>
class static_weighted_striped_square : public static_weighted_range<pos_t> {
public:
    using point_t = static_weighted_range<pos_t>::point_t;

protected:
    pos_t num_points = 0;
    pos_t seg_size = 0;
    std::vector<point_t> points;

public:
    static_weighted_striped_square() = default;

    static_weighted_striped_square(
        std::vector<point_t>& points,
        [[maybe_unused]] pos_t pos_max, [[maybe_unused]] uint16_t p = 1,
        pos_t seg_size = 128)
        : num_points(points.size())
        , seg_size(seg_size)
    {
        // this can also be done in O(n log seg_size) time, but requires more space
        ips4o::parallel::sort(points.begin(), points.end(),
            [&](const point_t& p1, const point_t& p2) {
                pos_t x1_s = p1.x / seg_size;
                pos_t x2_s = p2.x / seg_size;
                return x1_s < x2_s || (x1_s == x2_s && p1.y < p2.y);
            });

        this->points = std::move(points);
    }

    inline pos_t size() const
    {
        return num_points;
    }

    std::tuple<point_t, bool> lighter_point_in_range(
        pos_t weight,
        pos_t x1, pos_t x2,
        pos_t y1, pos_t y2) const
    {
        pos_t x1_s = x1 / seg_size;
        pos_t x2_s = x2 / seg_size;

        for (pos_t x = x1_s; x <= x2_s; x++) {
            pos_t b = x * seg_size;
            pos_t t = std::min<pos_t>(num_points, b + seg_size);
            pos_t m;

            while (b != t) {
                m = b + (t - b) / 2;
                const point_t& p = points[m];

                if (p.y < y1) {
                    b = m + 1;
                } else if (p.y > y2) {
                    t = m;
                } else
                    break;
            }

            if (b != t) {
                for (pos_t i = m; i < t; i++) {
                    const point_t& p = points[i];
                    if (p.y > y2)
                        break;

                    if (p.weight < weight && x1 <= p.x && p.x <= x2) {
                        return { p, true };
                    }
                }

                int64_t b_ = b;

                for (int64_t i = (int64_t)(m) - 1; i >= b_; i--) {
                    const point_t& p = points[i];
                    if (p.y < y1)
                        break;

                    if (p.weight < weight && x1 <= p.x && p.x <= x2) {
                        return { p, true };
                    }
                }
            }
        }

        return { { 0, 0, 0 }, false };
    }

    inline uint64_t size_in_bytes() const
    {
        return sizeof(this) + points.size() * sizeof(point_t);
    }

    static constexpr std::string name()
    {
        return "static weighted striped square";
    }
};

template <typename sidx_t = uint32_t>
using decomposed_static_weighted_striped_square =
    decomposed_range<static_weighted_striped_square, sidx_t>;