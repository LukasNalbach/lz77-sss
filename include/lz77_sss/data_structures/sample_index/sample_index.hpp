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

#include <cstdint>
#include <cstring>
#include <vector>

#include <ds/lce_naive_wordwise_xor.hpp>
#include <text/direct_text.hpp>
#include <lz77_sss/data_structures/rabin_karp_substring.hpp>
#include <lz77_sss/misc/log.hpp>
#include <lz77_sss/misc/search.hpp>
#include <lz77_sss/misc/utils.hpp>
#include <lz77_sss/data_structures/interval_hash_set.hpp>

enum class interval_samples {
    use,
    skip
};

template <typename text_t = lce::text::direct_text<char>,
    typename lce_r_t = lce::ds::lce_naive_wordwise_xor<char>,
    typename array_t = std::vector<uint40_t>>
class sample_index {
public:
    template <typename pos_t>
    struct interval {
        pos_t b;
        pos_t e;

        bool empty() const { return uint64_t(b) > uint64_t(e); }
    };

    using interval_t = interval<uint40_t>;

    struct query_ctx_t {
        uint64_t b;
        uint64_t e;
        uint64_t lce_b;
        uint64_t lce_e;

        inline uint64_t match_length() const { return std::min<uint64_t>(lce_b, lce_e); }

        inline interval_t interval() const
        {
            return { b, e };
        }
    };

protected:
    static constexpr uint64_t no_occ = largest_value<uint40_t>();
    static constexpr bool byte_text = text_t::is_byte_text;
    static constexpr uint64_t first_hashed_len_idx = byte_text ? 2 : 0;

    using rks_t = rabin_karp_substring<byte_text ? 31 : 61, text_t>;

    uint64_t n = 0;
    uint64_t s = 0;
    uint64_t max_patt_len_left = 0;
    uint64_t byte_size = 0;

    text_t T;
    array_view_t<array_t> S;
    const lce_r_t* LCE_R = nullptr;
    rks_t RKS;

    bit_aligned_vector PA_S;
    bit_aligned_vector SA_S;

    interval_t XIV_S_1[1 << 8];
    interval_t XIV_S_2[2][1 << 16];
    std::vector<interval_hash_set> XIV_S[2];

    std::vector<uint64_t> smpl_patt_lens[2];

    template <direction dir>
    inline uint64_t XA_S(uint64_t i) const
    {
        if constexpr (dir == LEFT) {
            return PA_S[i];
        } else {
            return SA_S[i];
        }
    }

    template <direction dir>
    inline bool is_pos_in_T(uint64_t p, uint64_t offs) const
    {
        if constexpr (dir == LEFT) {
            return p >= offs;
        } else {
            return p + offs < n;
        }
    }

    template <direction dir>
    inline uint16_t pair_at(uint64_t p) const
    {
        if constexpr (dir == LEFT) {
            return uint16_t(T[p - 1]) | (uint16_t(T[p]) << 8);
        } else {
            return uint16_t(T[p]) | (uint16_t(T[p + 1]) << 8);
        }
    }

    inline const interval_t& xiv_s_1(uint64_t pos_char) const { return XIV_S_1[T[pos_char]]; }

    template <direction dir>
    inline const interval_t& xiv_s_2(uint64_t pos_patt) const { return XIV_S_2[dir][pair_at<dir>(pos_patt)]; }

    inline bool occurs1(uint64_t pos_char) const { return xiv_s_1(pos_char).b != no_occ; }

    template <direction dir>
    inline bool occurs2(uint64_t pos_patt) const
    {
        if constexpr (dir == LEFT) {
            if (pos_patt < 1) [[unlikely]] return false;
        } else {
            if (pos_patt >= n - 1) [[unlikely]] return false;
        }

        return xiv_s_2<dir>(pos_patt).b != no_occ;
    }

    template <direction dir>
    inline uint64_t lce(
        uint64_t i, uint64_t j,
        uint64_t max_lce = std::numeric_limits<uint64_t>::max()) const
    {
        if constexpr (dir == LEFT) {
            return T.lce_left(i, j, max_lce);
        } else {
            return LCE_R->lce(i, j);
        }
    }

    template <direction dir>
    inline uint64_t lce_offs(
        uint64_t i, uint64_t j, uint64_t known_lce,
        uint64_t max_lce = std::numeric_limits<uint64_t>::max()) const
    {
        if constexpr (dir == LEFT) {
            if (!is_pos_in_T<dir>(std::min<uint64_t>(i, j), known_lce)) [[unlikely]] return known_lce;
            return known_lce + T.lce_left(i - known_lce, j - known_lce, max_lce - known_lce);
        } else {
            if (!is_pos_in_T<dir>(std::max<uint64_t>(i, j), known_lce)) [[unlikely]] return known_lce;
            return known_lce + LCE_R->lce(known_lce + i, known_lce + j);
        }
    }

    template <direction dir>
    inline bool cmp_lex(uint64_t i, uint64_t j, uint64_t common_len) const
    {
        if (i == j) [[unlikely]] return false;

        if constexpr (dir == LEFT) {
            if (common_len > std::min<uint64_t>(i, j)) [[unlikely]] return i < j;
            return T[i - common_len] < T[j - common_len];
        } else {
            if (std::max<uint64_t>(i, j) + common_len == n) [[unlikely]] return i > j;
            return T[i + common_len] < T[j + common_len];
        }
    }

    template <direction dir>
    inline bool cmp_sample_lex(uint64_t i, uint64_t j) const;

    void build_xiv_s_1_2(uint16_t p, bool log);

    template <direction dir>
    void build_interval_hash_sets(uint64_t max_smpl_len, uint16_t p, bool log);

public:
    sample_index() = default;

    void build(
        const text_t& T,
        uint64_t n,
        const array_t& S,
        const lce_r_t& LCE_R,
        interval_samples mode = interval_samples::use,
        uint64_t rks_sample_rate = 32,
        uint16_t p = 1,
        bool log = false,
        uint64_t max_patt_len_left = std::numeric_limits<uint64_t>::max(),
        uint64_t max_smpl_len_right = std::numeric_limits<uint64_t>::max());

    inline const rks_t& rks() const { return RKS; }

    inline uint64_t size_in_bytes() const { return byte_size; }

    inline uint64_t num_samples() const { return s; }

    inline uint64_t sample(uint64_t i) const { return S[i]; }

    inline uint64_t pa(uint64_t i) const { return PA_S[i]; }

    inline uint64_t sa(uint64_t i) const { return SA_S[i]; }

    inline query_ctx_t query() const
    {
        return {
            .b = 0,
            .e = s - 1,
            .lce_b = 0,
            .lce_e = 0
        };
    }

    template <direction dir>
    inline query_ctx_t query(const interval_t& xa_iv, uint64_t pos_patt, uint64_t len) const
    {
        return {
            .b = xa_iv.b,
            .e = xa_iv.e,
            .lce_b = lce_offs<dir>(pos_patt, S[XA_S<dir>(xa_iv.b)], len, max_patt_len_left),
            .lce_e = lce_offs<dir>(pos_patt, S[XA_S<dir>(xa_iv.e)], len, max_patt_len_left)
        };
    }

    inline query_ctx_t query_left(const interval_t& xa_iv, uint64_t pos_patt, uint64_t len) const
    {
        return query<LEFT>(xa_iv, pos_patt, len);
    }

    inline query_ctx_t query_right(const interval_t& xa_iv, uint64_t pos_patt, uint64_t len) const
    {
        return query<RIGHT>(xa_iv, pos_patt, len);
    }

    template <direction dir>
    inline uint64_t num_sampled_pattern_lengths() const { return smpl_patt_lens[dir].size(); }

    inline uint64_t num_sampled_pattern_lengths_left() const { return smpl_patt_lens[LEFT].size(); }

    inline uint64_t num_sampled_pattern_lengths_right() const
    {
        return smpl_patt_lens[RIGHT].size();
    }

    template <direction dir>
    inline const std::vector<uint64_t>& sampled_pattern_lengths() const
    {
        return smpl_patt_lens[dir];
    }

    inline const std::vector<uint64_t>& sampled_pattern_lengths_left() const
    {
        return smpl_patt_lens[LEFT];
    }

    inline const std::vector<uint64_t>& sampled_pattern_lengths_right() const
    {
        return smpl_patt_lens[RIGHT];
    }

    template <direction dir>
    inline std::pair<interval_t, bool> xa_interval(
        uint64_t patt_len_idx, uint64_t pos_patt,
        std::size_t hash = std::numeric_limits<std::size_t>::max()) const;

    inline std::pair<interval_t, bool> pa_interval(
        uint64_t patt_len_idx, uint64_t pos_patt,
        std::size_t hash = std::numeric_limits<std::size_t>::max()) const
    {
        return xa_interval<LEFT>(patt_len_idx, pos_patt, hash);
    }

    inline std::pair<interval_t, bool> sa_interval(
        uint64_t patt_len_idx, uint64_t pos_patt,
        std::size_t hash = std::numeric_limits<std::size_t>::max()) const
    {
        return xa_interval<RIGHT>(patt_len_idx, pos_patt, hash);
    }

    template <direction dir>
    inline bool extend(const query_ctx_t& qc_old, query_ctx_t& qc_new, uint64_t pos_patt, uint64_t len, interval_samples mode = interval_samples::use) const;

    bool extend_left(const query_ctx_t& qc_old, query_ctx_t& qc_new, uint64_t pos_patt, uint64_t len, interval_samples mode = interval_samples::use) const
    {
        return extend<LEFT>(qc_old, qc_new, pos_patt, len, mode);
    }

    inline bool extend_right(const query_ctx_t& qc_old, query_ctx_t& qc_new, uint64_t pos_patt, uint64_t len, interval_samples mode = interval_samples::use) const
    {
        return extend<RIGHT>(qc_old, qc_new, pos_patt, len, mode);
    }

    inline bool extend_left(query_ctx_t& qc, uint64_t pos_patt, uint64_t len, interval_samples mode = interval_samples::use) const
    {
        return extend<LEFT>(qc, qc, pos_patt, len, mode);
    }

    inline bool extend_right(query_ctx_t& qc, uint64_t pos_patt, uint64_t len, interval_samples mode = interval_samples::use) const
    {
        return extend<RIGHT>(qc, qc, pos_patt, len, mode);
    }

    template <direction dir>
    inline bool extend(query_ctx_t& qc, uint64_t pos_patt, uint64_t len, interval_samples mode = interval_samples::use) const
    {
        return extend<dir>(qc, qc, pos_patt, len, mode);
    }

    template <direction dir>
    query_ctx_t interpolate(const query_ctx_t& qc_short, const query_ctx_t& qc_long, uint64_t pos_patt, uint64_t len) const;

    inline query_ctx_t interpolate_left(const query_ctx_t& qc_short, const query_ctx_t& qc_long, uint64_t pos_patt, uint64_t len) const
    {
        return interpolate<LEFT>(qc_short, qc_long, pos_patt, len);
    };

    inline query_ctx_t interpolate_right(const query_ctx_t& qc_short, const query_ctx_t& qc_long, uint64_t pos_patt, uint64_t len) const
    {
        return interpolate<RIGHT>(qc_short, qc_long, pos_patt, len);
    };

    template <direction dir>
    inline void locate(const query_ctx_t& qc, std::vector<uint64_t>& Occ) const
    {
        Occ.reserve(Occ.size() + qc.e - qc.b + 1);

        for (uint64_t i = qc.b; i <= qc.e; i++) {
            Occ.emplace_back(S[XA_S<dir>(i)]);
        }
    }

    template <direction dir>
    inline std::vector<uint64_t> locate(const query_ctx_t& qc) const
    {
        std::vector<uint64_t> Occ;
        locate<dir>(qc, Occ);
        return Occ;
    }
};

#include "construction.cpp"
#include "queries.cpp"