
#include <algorithm>
#include <limits>
#include <vector>

#include <lz77_sss/data_structures/decomposed_range.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_range.hpp>
#include <lz77_sss/misc/utils.hpp>

template <typename pos_t = uint32_t>
class static_weighted_kd_tree : public static_weighted_range<pos_t> {
public:
    using point_t = static_weighted_range<pos_t>::point_t;

    static_weighted_kd_tree() = default;

    static_weighted_kd_tree(
        std::vector<point_t>& points,
        pos_t pos_max, uint16_t p = 1)
    {
        nodes.reserve(points.size());
        build<true>(points, 0, points.size());
        points.clear();
        points.shrink_to_fit();
    }

    std::tuple<point_t, bool> lighter_point_in_range(
        pos_t weight, pos_t x1, pos_t x2, pos_t y1, pos_t y2) const override
    {
        return lighter_point_in_range<true>(weight, 0, x1, x2, y1, y2);
    }

    inline pos_t size() const override
    {
        return nodes.size();
    }

    inline uint64_t size_in_bytes() const override
    {
        return sizeof(this) + nodes.size() * sizeof(node_t);
    }

protected:
    struct __attribute__((packed)) node_t {
        point_t point;
        pos_t left_child;
        pos_t right_child;
        pos_t min_weight;
    };

    std::vector<node_t> nodes;

    template <bool vertical>
    pos_t build(std::vector<point_t>& points, pos_t beg, pos_t end)
    {
        if (beg >= end) {
            return points.size();
        }

        pos_t mid = beg + (end - beg) / 2;

        std::nth_element(points.begin() + beg, points.begin() + mid, points.begin() + end,
            [](const point_t& a, const point_t& b) {
                return vertical ? a.x < b.x : a.y < b.y;
            });

        pos_t node_idx = nodes.size();

        nodes.emplace_back(node_t {
            .point = points[mid],
            .min_weight = points[mid].weight });

        node_t& node = nodes[node_idx];
        node.left_child = build<!vertical>(points, beg, mid);
        node.right_child = build<!vertical>(points, mid + 1, end);

        if (node.left_child != points.size()) {
            node.min_weight = std::min(
                node.min_weight,
                nodes[node.left_child].min_weight);
        }

        if (node.right_child != points.size()) {
            node.min_weight = std::min(
                node.min_weight,
                nodes[node.right_child].min_weight);
        }

        return node_idx;
    }

public:
    template <bool vertical>
    inline std::tuple<point_t, bool> lighter_point_in_range(
        pos_t weight, pos_t node_idx,
        pos_t x1, pos_t x2, pos_t y1, pos_t y2) const
    {
        if (node_idx == nodes.size()) {
            return { { 0, 0, 0 }, false };
        }

        const node_t& node = nodes[node_idx];
        const point_t& point = node.point;

        if (point.weight < weight && x1 <= point.x && point.x <= x2 && y1 <= point.y && point.y <= y2) {
            return { point, true };
        }

        if (((vertical && x1 <= point.x) || (!vertical && y1 <= point.y)) &&
            node.left_child != nodes.size() && nodes[node.left_child].min_weight < weight
        ) {
            auto [point_left, lighter_left] = lighter_point_in_range<!vertical>(
                weight, node.left_child, x1, x2, y1, y2);

            if (lighter_left) {
                return { point_left, true };
            }
        }

        if (((vertical && point.x <= x2) || (!vertical && point.y <= y2)) &&
            node.right_child != nodes.size() && nodes[node.right_child].min_weight < weight
        ) {
            auto [point_right, lighter_right] = lighter_point_in_range<!vertical>(
                weight, node.right_child, x1, x2, y1, y2);

            if (lighter_right) {
                return { point_right, true };
            }
        }

        return { { 0, 0, 0 }, false };
    }

    static constexpr std::string name()
    {
        return "static weighted k-d tree";
    }
};

template <typename sidx_t = uint32_t>
using decomposed_static_weighted_kd_tree =
    decomposed_range<static_weighted_kd_tree, sidx_t>;