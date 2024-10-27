#pragma once

#include <cstdint>

template <typename pos_t = uint32_t>
class static_weighted_range {
public:
    static constexpr bool is_decomposed() { return false; }
    static constexpr bool is_static() { return true; }
    static constexpr bool is_dynamic() { return false; }

    struct __attribute__((packed)) point_t {
        pos_t x, y;
        pos_t weight;
    };

    virtual pos_t size() const { }
    virtual uint64_t size_in_bytes() const { }

    virtual std::tuple<point_t, bool> lighter_point_in_range(
        pos_t weight, pos_t x1, pos_t x2, pos_t y1, pos_t y2) const { }
};