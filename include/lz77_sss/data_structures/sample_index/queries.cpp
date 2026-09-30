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

#include <lz77_sss/data_structures/sample_index/sample_index.hpp>

template <typename text_t, typename lce_r_t, typename array_t>
template <direction dir>
inline std::pair<typename sample_index<text_t, lce_r_t, array_t>::interval_t, bool>
sample_index<text_t, lce_r_t, array_t>::xa_interval(
    uint64_t patt_len_idx, uint64_t pos_patt, std::size_t hash) const
{
    if constexpr (byte_text) {
        if (patt_len_idx == 0) {
            if (!occurs1(pos_patt)) [[unlikely]] {
                return { { 0, 0 }, false };
            } else {
                return { xiv_s_1(pos_patt), true };
            }
        } else if (patt_len_idx == 1) {
            if (!occurs2<dir>(pos_patt)) [[unlikely]] {
                return { { 0, 0 }, false };
            } else {
                return { xiv_s_2<dir>(pos_patt), true };
            }
        }
    }

    const uint64_t patt_len = sampled_pattern_lengths<dir>()[patt_len_idx];

    if constexpr (!byte_text) {
        if (!is_pos_in_T<dir>(pos_patt, patt_len - 1)) [[unlikely]] return { { 0, 0 }, false };
    }

    if (hash == std::numeric_limits<std::size_t>::max()) {
        hash = RKS.template substr_fp<dir>(pos_patt, patt_len);
    }

    uint64_t b;
    uint64_t e;

    const bool found = XIV_S[dir][patt_len_idx].find(hash, [&](uint64_t beg) {
        const uint64_t pos = sample(XA_S<dir>(beg));
        return pos == pos_patt || lce<dir>(pos, pos_patt, patt_len) >= patt_len;
    }, b, e);

    if (!found) {
        return { { 0, 0 }, false };
    }

    return { interval_t { .b = uint40_t(b), .e = uint40_t(e) }, true };
}

template <typename text_t, typename lce_r_t, typename array_t>
template <direction dir>
bool sample_index<text_t, lce_r_t, array_t>::extend(
    const query_ctx_t& qc_old, query_ctx_t& qc_new,
    uint64_t pos_patt, uint64_t len, interval_samples mode) const
{
    if (qc_old.match_length() >= len) {
        if (&qc_old != &qc_new) {
            qc_new = qc_old;
        }

        return true;
    }

    const bool use_interval_samples = mode == interval_samples::use && num_sampled_pattern_lengths<dir>() > 0;

    uint64_t b;
    uint64_t e;
    uint64_t lce_b;
    uint64_t lce_e;

    uint64_t patt_len_idx = 0;
    uint64_t smpl_len;
    std::size_t fp_smpl = std::numeric_limits<std::size_t>::max();

    if (use_interval_samples) {
        patt_len_idx = bin_search_max_leq<uint64_t, uint64_t>(
            len, 0, num_sampled_pattern_lengths<dir>() - 1, [&](uint64_t x) {
                return sampled_pattern_lengths<dir>()[x];
            });

        smpl_len = sampled_pattern_lengths<dir>()[patt_len_idx];
    }

    if (use_interval_samples && std::min<uint64_t>(qc_old.lce_b, qc_old.lce_e) < smpl_len) {
        if (patt_len_idx >= first_hashed_len_idx) {
            fp_smpl = RKS.template substr_fp<dir>(pos_patt, smpl_len);
        }

        auto [iv, found] = xa_interval<dir>(patt_len_idx, pos_patt, fp_smpl);

        if (!found) {
            return false;
        }

        b = iv.b;
        e = iv.e;
        lce_b = smpl_len;
        lce_e = smpl_len;
    } else {
        b = qc_old.b;
        e = qc_old.e;
        lce_b = qc_old.lce_b;
        lce_e = qc_old.lce_e;
    }

    if (lce_b < len) [[likely]]
        lce_b = lce_offs<dir>(pos_patt, S[XA_S<dir>(b)], lce_b, len);

    if (lce_e < len) [[likely]] {
        if (e == b) {
            lce_e = lce_b;
        } else {
            lce_e = lce_offs<dir>(pos_patt, S[XA_S<dir>(e)], lce_e, len);
        }
    }

    if (use_interval_samples && std::min<uint64_t>(lce_b, lce_e) < len &&
        patt_len_idx + 1 < num_sampled_pattern_lengths<dir>() &&
        is_pos_in_T<dir>(pos_patt, sampled_pattern_lengths<dir>()[patt_len_idx + 1] - 1)
    ) {
        uint64_t nxt_smpl_len = sampled_pattern_lengths<dir>()[patt_len_idx + 1];
        uint64_t len_diff = nxt_smpl_len - smpl_len;
        std::size_t fp_nxt_smpl = std::numeric_limits<std::size_t>::max();

        if (patt_len_idx + 1 >= first_hashed_len_idx) {
            if (fp_smpl == std::numeric_limits<std::size_t>::max()) {
                fp_nxt_smpl = RKS.template substr_fp<dir>(pos_patt, nxt_smpl_len);
            } else if constexpr (dir == LEFT) {
                fp_nxt_smpl = RKS.concat(
                    RKS.template substr_fp<LEFT>(pos_patt - smpl_len, len_diff),
                    fp_smpl, smpl_len);
            } else {
                fp_nxt_smpl = RKS.concat(
                    fp_smpl, RKS.template substr_fp<RIGHT>(pos_patt + smpl_len, len_diff),
                    len_diff);
            }
        }

        auto [iv2, found_nxt] = xa_interval<dir>(
            patt_len_idx + 1, pos_patt, fp_nxt_smpl);

        if (found_nxt) {
            qc_new = interpolate<dir>(
                { b, e, lce_b, lce_e },
                { iv2.b, iv2.e, nxt_smpl_len, nxt_smpl_len },
                pos_patt, len);

            return true;
        }
    }

    uint64_t e_min = b;
    uint64_t e_max = e;
    uint64_t lce_e_min = lce_b;
    uint64_t lce_e_max = lce_e;

    if (lce_b < len) {
        uint64_t lo = b;
        uint64_t hi = e;
        uint64_t lce_lo = lce_b;
        uint64_t lce_hi = lce_e;

        uint64_t mid;
        uint64_t pos_mid, lce_mid;

        while (hi - lo > 1) {
            mid = lo + (hi - lo) / 2;
            pos_mid = S[XA_S<dir>(mid)];

            lce_mid = lce_offs<dir>(
                pos_patt, pos_mid,
                std::min<uint64_t>(lce_lo, lce_hi),
                len);

            if (lce_mid >= len) {
                hi = mid;
                lce_hi = lce_mid;

                if (mid > e_min) {
                    e_min = mid;
                    lce_e_min = lce_mid;
                }
            } else if (cmp_lex<dir>(pos_mid, pos_patt, lce_mid)) {
                lo = mid;
                lce_lo = lce_mid;

                if (mid > e_min) {
                    e_min = mid;
                    lce_e_min = lce_mid;
                }
            } else {
                hi = mid;
                lce_hi = lce_mid;

                if (mid < e_max) {
                    e_max = mid;
                    lce_e_max = lce_mid;
                }
            }
        }

        if (lce_lo < len) {
            if (lce_hi < len) {
                return false;
            }

            qc_new.b = hi;
            qc_new.lce_b = lce_hi;
        } else {
            qc_new.b = lo;
            qc_new.lce_b = lce_lo;
        }
    } else {
        qc_new.b = b;
        qc_new.lce_b = lce_b;
    }

    if (lce_e < len) {
        uint64_t lo = e_min;
        uint64_t hi = e_max;
        uint64_t lce_lo = lce_e_min;
        uint64_t lce_hi = lce_e_max;

        uint64_t mid;
        uint64_t pos_mid, lce_mid;

        while (hi - lo > 1) {
            mid = lo + (hi - lo) / 2;
            pos_mid = S[XA_S<dir>(mid)];

            lce_mid = lce_offs<dir>(
                pos_patt, pos_mid,
                std::min<uint64_t>(lce_lo, lce_hi),
                len);

            if (lce_mid >= len) {
                lo = mid;
                lce_lo = lce_mid;
            } else {
                hi = mid;
                lce_hi = lce_mid;
            }
        }

        if (lce_hi >= len) {
            qc_new.e = hi;
            qc_new.lce_e = lce_hi;
        } else {
            qc_new.e = lo;
            qc_new.lce_e = lce_lo;
        }
    } else {
        qc_new.e = e;
        qc_new.lce_e = lce_e;
    }

    return true;
}

template <typename text_t, typename lce_r_t, typename array_t>
template <direction dir>
sample_index<text_t, lce_r_t, array_t>::query_ctx_t
sample_index<text_t, lce_r_t, array_t>::interpolate(
    const query_ctx_t& qc_short,
    const query_ctx_t& qc_long,
    uint64_t pos_patt, uint64_t len) const
{
    if (qc_short.match_length() >= len) {
        return qc_short;
    }

    uint64_t lo = qc_short.b;
    uint64_t hi = qc_long.b;
    uint64_t mid;

    uint64_t lce_lo = qc_short.lce_b;
    uint64_t lce_hi = qc_long.lce_b;
    uint64_t lce_mid;

    query_ctx_t qc_ret;
    uint64_t pos_mid;

    lce_lo = lce_offs<dir>(
        pos_patt, S[XA_S<dir>(lo)],
        lce_lo,
        len);

    while (hi - lo > 1) {
        mid = lo + (hi - lo) / 2;
        pos_mid = S[XA_S<dir>(mid)];

        lce_mid = lce_offs<dir>(
            pos_patt, pos_mid,
            std::min<uint64_t>(lce_lo, lce_hi),
            len);

        if (lce_mid < len) {
            lo = mid;
            lce_lo = lce_mid;
        } else {
            hi = mid;
            lce_hi = lce_mid;
        }
    }

    if (lce_lo < len) {
        qc_ret.b = hi;
        qc_ret.lce_b = lce_hi;
    } else {
        qc_ret.b = lo;
        qc_ret.lce_b = lce_lo;
    }

    lo = qc_long.e;
    hi = qc_short.e;

    lce_lo = qc_long.lce_e;
    lce_hi = qc_short.lce_e;

    lce_hi = lce_offs<dir>(
        pos_patt, S[XA_S<dir>(hi)],
        lce_hi, len);

    while (hi - lo > 1) {
        mid = lo + (hi - lo) / 2;
        pos_mid = S[XA_S<dir>(mid)];

        lce_mid = lce_offs<dir>(
            pos_patt, pos_mid,
            std::min<uint64_t>(lce_lo, lce_hi),
            len);

        if (lce_mid < len) {
            hi = mid;
            lce_hi = lce_mid;
        } else {
            lo = mid;
            lce_lo = lce_mid;
        }
    }

    if (lce_hi < len) {
        qc_ret.e = lo;
        qc_ret.lce_e = lce_lo;
    } else {
        qc_ret.e = hi;
        qc_ret.lce_e = lce_hi;
    }

    return qc_ret;
}