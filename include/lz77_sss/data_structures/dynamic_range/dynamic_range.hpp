#pragma once

#include <cstdint>

template<typename pos_t = uint32_t>
class dynamic_range {
    public:

    struct point_t {
        pos_t x, y;
    };

    virtual pos_t size() const {}
    virtual void insert(point_t p) {}
    virtual uint64_t size_in_bytes() const {}

    virtual std::tuple<point_t, bool> point_in_range(
        pos_t x1, pos_t x2, pos_t y1, pos_t y2
    ) const {}
};