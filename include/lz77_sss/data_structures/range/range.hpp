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

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <lz77_sss/data_structures/bit_aligned_interleaved_vectors.hpp>
#include <lz77_sss/data_structures/range/semi_dynamic_square_grid.hpp>
#include <lz77_sss/data_structures/range/static_weighted_kd_tree.hpp>
#include <lz77_sss/data_structures/range/static_weighted_square_grid.hpp>
#include <lz77_sss/data_structures/range/range_point.hpp>
#include <lz77_sss/misc/utils.hpp>

class range_ds {
public:
    using point_t = range_point_t;
    using result_t = std::tuple<point_t, bool>;

    virtual ~range_ds() = default;

    virtual bool is_decomposed() const = 0;
    virtual bool is_static() const = 0;
    virtual bool is_dynamic() const = 0;
    virtual std::string name() const = 0;
    virtual uint64_t size() const = 0;
    virtual uint64_t size_in_bytes() const = 0;
    virtual void insert(uint64_t c, point_t point) = 0;

    virtual result_t lighter_point_in_range(uint64_t c, uint64_t weight,
        uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const = 0;

    virtual result_t point_in_range(uint64_t c,
        uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const = 0;
};

template <typename impl_t>
class range_ds_impl final : public range_ds {
    impl_t impl;

public:
    range_ds_impl(std::vector<point_t>& points, uint64_t pos_max_excl, uint16_t p)
        : impl(points, pos_max_excl, p) { }

    bool is_decomposed() const override { return false; }
    bool is_static() const override { return impl_t::is_static(); }
    bool is_dynamic() const override { return impl_t::is_dynamic(); }
    std::string name() const override { return impl_t::name(); }
    uint64_t size() const override { return impl.size(); }
    uint64_t size_in_bytes() const override { return impl.size_in_bytes(); }

    void insert(uint64_t, point_t point) override { if constexpr (impl_t::is_dynamic()) impl.insert(point); }

    result_t lighter_point_in_range(uint64_t, uint64_t weight,
        uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const override
    {
        if constexpr (impl_t::is_static())
            return impl.lighter_point_in_range(weight, x1, x2, y1, y2);
        else
            return { point_t { }, false };
    }

    result_t point_in_range(uint64_t,
        uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2) const override
    {
        if constexpr (impl_t::is_dynamic())
            return impl.point_in_range(x1, x2, y1, y2);
        else
            return { point_t { }, false };
    }
};

template <typename impl_t>
class decomposed_range;

template <typename impl_t>
class grouped_range;

enum class range_ds_type {
    swsg,
    swkdt,
    sdsg
};

template <typename fnc_t>
inline auto visit_range_ds_type(range_ds_type type, fnc_t fnc)
{
    if (type == range_ds_type::swkdt) return fnc(std::type_identity<static_weighted_kd_tree>());
    if (type == range_ds_type::sdsg) return fnc(std::type_identity<semi_dynamic_square_grid>());
    return fnc(std::type_identity<static_weighted_square_grid>());
}

struct range_ds_kind {
    range_ds_type type = range_ds_type::swsg;
    bool decomposed = false;

    bool is_static() const
    {
        return visit_range_ds_type(type, [](auto impl) { return decltype(impl)::type::is_static(); });
    }

    bool is_dynamic() const { return !is_static(); }

    std::string name() const
    {
        const std::string inner = visit_range_ds_type(type, [](auto impl) { return decltype(impl)::type::name(); });
        return decomposed ? "d-" + inner : inner;
    }
};

inline std::optional<range_ds_kind> range_ds_by_name(const std::string& name)
{
    for (range_ds_type type : { range_ds_type::swsg, range_ds_type::swkdt, range_ds_type::sdsg }) {
        const range_ds_kind kind { .type = type, .decomposed = name.starts_with("d-") };
        if (kind.name() == name) return kind;
    }

    return std::nullopt;
}

inline range_ds* make_range_ds(const range_ds_kind& kind,
    std::vector<range_point_t>& points, uint64_t pos_max_excl, uint16_t p)
{
    return visit_range_ds_type(kind.type, [&](auto impl) -> range_ds* {
        return new range_ds_impl<typename decltype(impl)::type>(points, pos_max_excl, p);
    });
}

template <typename text_t, typename array_t, typename points_t>
inline range_ds* make_range_ds(const range_ds_kind& kind,
    const text_t& T, const array_t& S,
    const points_t& points, uint16_t p)
{
    constexpr bool byte_chars = sizeof(std::remove_cvref_t<decltype(T[0])>) == 1;
    std::vector<std::conditional_t<byte_chars, uint8_t, uint32_t>> chr;
    no_init_resize(chr, S.size());

    #pragma omp parallel for num_threads(p) schedule(dynamic, 65536)
    for (uint64_t i = 0; i < S.size(); i++) {
        chr[i] = T[S[i]];
    }

    bit_aligned_interleaved_vectors<3> packed;
    const bit_aligned_interleaved_vectors<3>* P = &packed;

    if constexpr (std::is_same_v<points_t, bit_aligned_interleaved_vectors<3>>) {
        P = &points;
    } else {
        std::array<uint64_t, 3> max_values = { 0, 0, 0 };

        for (const range_point_t& point : points) {
            max_values[0] = std::max<uint64_t>(max_values[0], point.x);
            max_values[1] = std::max<uint64_t>(max_values[1], point.y);
            max_values[2] = std::max<uint64_t>(max_values[2], point.weight);
        }

        packed = bit_aligned_interleaved_vectors<3>(points.size(), max_values);

        for (uint64_t i = 0; i < points.size(); i++) {
            packed.set<0>(i, points[i].x);
            packed.set<1>(i, points[i].y);
            packed.set<2>(i, points[i].weight);
        }
    }

    return visit_range_ds_type(kind.type, [&](auto impl) -> range_ds* {
        if constexpr (byte_chars) {
            return new decomposed_range<typename decltype(impl)::type>(chr, *P, p);
        } else {
            return new grouped_range<typename decltype(impl)::type>(chr, *P, p);
        }
    });
}

#include <lz77_sss/data_structures/range/decomposed_range.hpp>
#include <lz77_sss/data_structures/range/grouped_range.hpp>
