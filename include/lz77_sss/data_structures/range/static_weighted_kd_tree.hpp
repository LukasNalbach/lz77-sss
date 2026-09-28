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

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include <lz77_sss/data_structures/range/range_point.hpp>
#include <lz77_sss/misc/parallel_nth_element.hpp>
#include <lz77_sss/misc/utils.hpp>

class static_weighted_kd_tree {
    enum class axis {
        x,
        y
    };

    static constexpr axis other(axis a) { return a == axis::x ? axis::y : axis::x; }

public:
    static constexpr bool is_decomposed() { return false; }
    static constexpr bool is_static() { return true; }
    static constexpr bool is_dynamic() { return false; }

    using point_t = range_point_t;

    static_weighted_kd_tree() = default;

    static_weighted_kd_tree(
        std::vector<point_t>& points,
        [[maybe_unused]] uint64_t pos_max_excl, uint16_t p = 1)
    {
        no_init_resize(nodes, points.size());

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < points.size(); i++) {
            nodes[i] = { };
        }

        std::vector<subtree_t> top_nodes;
        std::vector<subtree_t> subtrees;
        axis split_axis = axis::x;

        if (!points.empty()) {
            subtrees.emplace_back(0, points.size(), 0);
        }

        while (!subtrees.empty() && subtrees.size() < p) {
            std::vector<subtree_t> children;

            for (const auto& [beg, end, node_idx] : subtrees) {
                const uint64_t mid = beg + (end - beg) / 2;

                if (split_axis == axis::x) {
                    parallel_nth_element(points, beg, mid, end, p, coord<axis::x>);
                } else {
                    parallel_nth_element(points, beg, mid, end, p, coord<axis::y>);
                }

                nodes[node_idx] = {
                    .point = points[mid],
                    .min_weight = points[mid].weight };

                if (beg < mid) children.emplace_back(beg, mid, node_idx + 1);
                if (mid + 1 < end) children.emplace_back(mid + 1, end, node_idx + 1 + (end - beg) / 2);
            }

            top_nodes.insert(top_nodes.end(), subtrees.begin(), subtrees.end());
            subtrees = std::move(children);
            split_axis = other(split_axis);
        }

        #pragma omp parallel for schedule(dynamic, 1) num_threads(p)
        for (uint64_t i = 0; i < subtrees.size(); i++) {
            const auto [beg, end, node_idx] = subtrees[i];

            if (split_axis == axis::x) {
                build<axis::x>(points, beg, end, node_idx);
            } else {
                build<axis::y>(points, beg, end, node_idx);
            }
        }

        for (uint64_t i = top_nodes.size(); i > 0; i--) {
            const auto [beg, end, node_idx] = top_nodes[i - 1];
            update_min_weight(node_idx, end - beg);
        }

        points.clear();
        points.shrink_to_fit();
    }

    std::tuple<point_t, bool> lighter_point_in_range(
        uint64_t weight, uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const
    {
        return lighter_point_in_range<axis::x>(weight, 0, size(), x1, x2, y1, y2);
    }

    inline uint64_t size() const { return nodes.size(); }

    inline uint64_t size_in_bytes() const { return sizeof(*this) + nodes.size() * sizeof(node_t); }

protected:
    struct __attribute__((packed)) node_t {
        point_t point;
        uint40_t min_weight;
    };

    using subtree_t = std::tuple<uint64_t, uint64_t, uint64_t>;

    static constexpr uint64_t min_par_chunk_size = 1 << 14;

    std::vector<node_t> nodes;

    template <axis split_axis>
    static inline uint64_t coord(const point_t& point) { return split_axis == axis::x ? point.x : point.y; }

    template <axis split_axis>
    void build(std::vector<point_t>& points, uint64_t beg, uint64_t end, uint64_t node_idx)
    {
        if (beg >= end) return;
        uint64_t mid = beg + (end - beg) / 2;

        std::nth_element(points.begin() + beg, points.begin() + mid, points.begin() + end,
            [](const point_t& a, const point_t& b) {
                return split_axis == axis::x ? a.x < b.x : a.y < b.y;
            });

        nodes[node_idx] = {
            .point = points[mid],
            .min_weight = points[mid].weight };

        build<other(split_axis)>(points, beg, mid, node_idx + 1);
        build<other(split_axis)>(points, mid + 1, end, node_idx + 1 + (end - beg) / 2);
        update_min_weight(node_idx, end - beg);
    }

    void update_min_weight(uint64_t node_idx, uint64_t subtree_size)
    {
        node_t& node = nodes[node_idx];
        const uint64_t left_child = node_idx + 1;
        const uint64_t right_child = left_child + subtree_size / 2;
        const bool has_left_child = subtree_size > 1;
        const bool has_right_child = subtree_size > 2;

        if (has_left_child) {
            node.min_weight = std::min(
                node.min_weight,
                nodes[left_child].min_weight);
        }

        if (has_right_child) {
            node.min_weight = std::min(
                node.min_weight,
                nodes[right_child].min_weight);
        }
    }

public:
    template <axis split_axis>
    inline std::tuple<point_t, bool> lighter_point_in_range(
        uint64_t weight, uint64_t node_idx, uint64_t subtree_end,
        uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const
    {
        const node_t& node = nodes[node_idx];
        const point_t& point = node.point;
        const uint64_t left_child = node_idx + 1;
        const uint64_t right_child = left_child + (subtree_end - node_idx) / 2;
        const bool has_left_child = subtree_end - node_idx > 1;
        const bool has_right_child = subtree_end - node_idx > 2;

        if (point.weight < weight &&
            x1 <= point.x && point.x <= x2 &&
            y1 <= point.y && point.y <= y2
        ) {
            return { point, true };
        }

        if (has_left_child && nodes[left_child].min_weight < weight &&
            ((split_axis == axis::x && x1 <= point.x) || (split_axis == axis::y && y1 <= point.y))
        ) {
            auto [point_left, found] = lighter_point_in_range<other(split_axis)>(
                weight, left_child, right_child, x1, x2, y1, y2);

            if (found) {
                return { point_left, true };
            }
        }

        if (has_right_child && nodes[right_child].min_weight < weight &&
            ((split_axis == axis::x && point.x <= x2) || (split_axis == axis::y && point.y <= y2))
        ) {
            auto [point_right, found] = lighter_point_in_range<other(split_axis)>(
                weight, right_child, subtree_end, x1, x2, y1, y2);

            if (found) {
                return { point_right, true };
            }
        }

        return { { 0, 0, 0 }, false };
    }

    static constexpr std::string name() { return "swkdt"; }
};

