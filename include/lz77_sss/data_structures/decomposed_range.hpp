#pragma once

#include <type_traits>

#include <lz77_sss/misc/utils.hpp>
#include <lz77_sss/data_structures/dynamic_range/dynamic_range.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_range.hpp>

template <template <typename> typename range_ds_t, typename sidx_t>
using range_ds_base = std::conditional_t<
    std::is_base_of_v<static_weighted_range<sidx_t>, range_ds_t<sidx_t>>,
    static_weighted_range<sidx_t>, dynamic_range<sidx_t>>;

template <template <typename> typename range_ds_t, typename sidx_t = uint32_t, typename char_t = char>
class decomposed_range : public range_ds_base<range_ds_t, sidx_t> {
public:
    using point_t = range_ds_t<sidx_t>::point_t;
    static constexpr bool is_decomposed() { return true; }
    static constexpr bool is_static() { return range_ds_t<sidx_t>::is_static(); }
    static constexpr bool is_dynamic() { return range_ds_t<sidx_t>::is_dynamic(); }

protected:
    sidx_t num_points = 0;
    std::array<sidx_t, 257> C_S = { 0 };
    std::array<range_ds_t<sidx_t>, 256> R_c;

    inline void to_internal(uint8_t uc, point_t& p) const
    {
        sidx_t rnk_c = C_S[uc];
        p.x -= rnk_c;
        p.y -= rnk_c;
    }

    inline void to_internal(uint8_t uc,
        sidx_t& x1, sidx_t& x2,
        sidx_t& y1, sidx_t& y2) const
    {
        sidx_t rnk_c = C_S[uc];
        x1 -= rnk_c;
        x2 -= rnk_c;
        y1 -= rnk_c;
        y2 -= rnk_c;
    }

    inline void to_external(uint8_t uc, point_t& p) const
    {
        sidx_t rnk_c = C_S[uc];
        p.x += rnk_c;
        p.y += rnk_c;
    }

public:
    decomposed_range() = default;

    template <typename pos_t>
    decomposed_range(
        const char_t* T,
        const std::vector<pos_t>& S,
        const std::vector<point_t>& P,
        uint16_t p = 1)
    {
        std::array<std::vector<point_t>, 256> P_c;
        std::vector<std::vector<pos_t>> C_p(p,
            std::vector<pos_t>(256, 0));

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < S.size(); i++) {
            C_p[omp_get_thread_num()][char_to_uchar(T[S[i]])]++;
        }

        for (uint16_t c = 0; c < 256; c++) {
            for (uint16_t j = 0; j < p; j++) {
                C_S[c] += C_p[j][c];
            }
        }
        
        for (uint16_t c = 1; c < 256; c++) C_S[c] += C_S[c - 1];
        for (uint16_t c = 256; c > 0; c--) C_S[c] = C_S[c - 1];
        C_S[0] = 0;

        #pragma omp parallel for num_threads(p)
        for (uint16_t c = 0; c < 256; c++) {
            P_c[c].reserve(C_S[c + 1] - C_S[c]);
        }

        for (sidx_t i = 0; i < S.size(); i++) {
            uint8_t uc = char_to_uchar(T[S[i]]);
            point_t p = P[i];
            to_internal(uc, p);
            P_c[uc].emplace_back(p);
        }

        for (uint16_t c = 0; c < 256; c++) {
            sidx_t frq = C_S[c + 1] - C_S[c];
            if (frq == 0) continue;
            R_c[c] = range_ds_t<sidx_t>(P_c[c], frq, p);
            P_c[c].clear();
            P_c[c].shrink_to_fit();
        }

        if constexpr (is_static()) {
            num_points = S.size();
        }
    }

    inline sidx_t size() const override
    {
        return num_points;
    }

    inline sidx_t size(char c) const
    {
        return R_c[char_to_uchar(c)].size();
    }

    inline uint64_t size_in_bytes() const override
    {
        uint64_t bytes = sizeof(this);

        for (uint16_t c = 0; c < 256; c++) {
            bytes += R_c[c].size_in_bytes();
        }

        return bytes;
    }

    inline void insert(char c, point_t p)
        requires(is_dynamic())
    {
        uint8_t uc = char_to_uchar(c);
        to_internal(uc, p);
        R_c[uc].insert(p);
        num_points++;
    }

    inline std::tuple<point_t, bool> lighter_point_in_range(
        char c, sidx_t weight,
        sidx_t x1, sidx_t x2,
        sidx_t y1, sidx_t y2) const
        requires(is_static())
    {
        uint8_t uc = char_to_uchar(c);
        to_internal(uc, x1, x2, y1, y2);
        auto [p, res] = R_c[uc].lighter_point_in_range(
            weight, x1, x2, y1, y2);
        to_external(uc, p);
        return { p, res };
    }

    inline std::tuple<point_t, bool> point_in_range(
        char c,
        sidx_t x1, sidx_t x2,
        sidx_t y1, sidx_t y2) const
        requires(is_dynamic())
    {
        uint8_t uc = char_to_uchar(c);
        to_internal(uc, x1, x2, y1, y2);
        auto [p, res] = R_c[uc].point_in_range(x1, x2, y1, y2);
        to_external(uc, p);
        return { p, res };
    }

    static constexpr std::string name()
    {
        return "decomposed " + range_ds_t<sidx_t>::name();
    }
};