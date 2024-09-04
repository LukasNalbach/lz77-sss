#pragma once

template <typename pos_t, typename sidx_t, typename lce_r_t>
template <direction dir>
inline std::pair<typename sample_index<pos_t, sidx_t, lce_r_t>::sxa_interval_t, bool>
sample_index<pos_t, sidx_t, lce_r_t>::sxa_interval(
    pos_t pat_len_idx, pos_t pos_pat, std::size_t hash
) {
    if (pat_len_idx == 0) {
        if (!occurs1(pos_pat)) {
            return {{0, 0}, false};
        } else {
            return {sciv(pos_pat), true};
        }
    } else if (pat_len_idx == 1) {
        if (!occurs2<dir>(pos_pat)) {
            return {{0, 0}, false};
        } else {
            return {sxiv2<dir>(pos_pat), true};
        }
    }

    const auto& map = SXIVX<dir>()[pat_len_idx];

    if (hash == std::numeric_limits<std::size_t>::max()) {
        pos_t pat_len = sampled_pattern_lengths<dir>()[pat_len_idx];
        hash = rks.substring<dir>(pos_pat, pat_len);
    }

    pos_hash = pos_pat;
    auto it = map.find({s, s}, hash);

    if (it == map.end()) {
        return {{0, 0}, false};
    } else {
        return {*it, true};
    }
}

template <typename pos_t, typename sidx_t, typename lce_r_t>
template <direction dir>
bool sample_index<pos_t, sidx_t, lce_r_t>::extend(
    const query_context& qc_old, query_context& qc_new,
    pos_t pos_pat, pos_t len, bool use_samples
) {
    if (qc_old.match_length() >= len) {
        if (&qc_old != &qc_new) {
            qc_new = qc_old;
        }

        return true;
    }

    sidx_t b;
    sidx_t e;
    pos_t lce_b;
    pos_t lce_e;

    pos_t pat_len_idx;
    pos_t len_smpl;
    std::size_t fp_smpl = std::numeric_limits<std::size_t>::max();

    if (use_samples) {
        pat_len_idx = bin_search_max_leq<pos_t, pos_t>(
            len, 0, num_sampled_pattern_lengths<dir>() - 1, [&](pos_t x){
                return sampled_pattern_lengths<dir>()[x];
        });

        len_smpl = sampled_pattern_lengths<dir>()[pat_len_idx];
    }
    
    if (use_samples && std::min<pos_t>(qc_old.lce_b, qc_old.lce_e) < len_smpl) {
        if (pat_len_idx >= 2) {
            fp_smpl = rks.substring<dir>(pos_pat, len_smpl);
        }

        auto [iv, result] = sxa_interval<dir>(pat_len_idx, pos_pat, fp_smpl);

        if (!result) {
            return false;
        }

        b = iv.b;
        e = iv.e;
        lce_b = len_smpl;
        lce_e = len_smpl;
    } else {
        b = qc_old.b;
        e = qc_old.e;
        lce_b = qc_old.lce_b;
        lce_e = qc_old.lce_e;
    }

    if (lce_b < len) lce_b = lce_offs<dir>(
        pos_pat, S[SXA<dir>(b)], lce_b, len);
    
    if (lce_e < len) {
        if (e == b) {
            lce_e = lce_b;
        } else {
            lce_e = lce_offs<dir>(
                pos_pat, S[SXA<dir>(e)], lce_e, len);
        }
    }

    if (use_samples && std::min<pos_t>(lce_b, lce_e) < len &&
        pat_len_idx + 1 < num_sampled_pattern_lengths<dir>() &&
        is_pos_in_T<dir>(pos_pat, sampled_pattern_lengths<dir>()[pat_len_idx + 1] - 1)
    ) {
        pos_t len_nxt_smpl = sampled_pattern_lengths<dir>()[pat_len_idx + 1];
        pos_t len_diff = len_nxt_smpl - len_smpl;
        std::size_t fp_nxt_smpl;

        if (pat_len_idx >= 1) {
            if (fp_smpl == std::numeric_limits<std::size_t>::max()) {
                fp_nxt_smpl = rks.substring<dir>(pos_pat, len_nxt_smpl);
            } else {
                if constexpr (dir == LEFT) {
                    fp_nxt_smpl = rks.concat(
                        rks.substring<LEFT>(pos_pat - len_smpl, len_diff),
                        fp_smpl, len_smpl);
                } else {
                    fp_nxt_smpl = rks.concat(
                        fp_smpl, rks.substring<RIGHT>(
                            pos_pat + len_smpl, len_diff),
                        len_diff);
                }
            }
        }
        
        auto [iv2, result2] = sxa_interval<dir>(
            pat_len_idx + 1, pos_pat, fp_nxt_smpl);
        
        if (result2) {
            qc_new = interpolate<dir>(
                {b, e, lce_b, lce_e},
                {iv2.b, iv2.e, len_nxt_smpl, len_nxt_smpl},
                pos_pat, len
            );

            return true;
        }
    }
    
    sidx_t e_min = b;
    sidx_t e_max = e;
    pos_t lce_e_min = lce_b;
    pos_t lce_e_max = lce_e;

    if (lce_b < len) {
        sidx_t l = b;
        sidx_t r = e;
        pos_t lce_l = lce_b;
        pos_t lce_r = lce_e;

        sidx_t m;
        pos_t pos_m, lce_m;

        while (r - l > 1) {
            m = l + (r - l) / 2;
            pos_m = S[SXA<dir>(m)];

            lce_m = lce_offs<dir>(
                pos_pat, pos_m,
                std::min<pos_t>(lce_l, lce_r),
                len
            );

            if (lce_m >= len) {
                r = m;
                lce_r = lce_m;

                if (m > e_min) {
                    e_min = m;
                    lce_e_min = lce_m;
                }
            } else if (cmp_lex<dir>(
                pos_m, pos_pat, lce_m
            )) {
                l = m;
                lce_l = lce_m;
            } else {
                r = m;
                lce_r = lce_m;

                if (m < e_max) {
                    e_max = m;
                    lce_e_max = lce_m;
                }
            }
        }

        if (lce_l < len) {
            if (lce_r < len) {
                return false;
            }

            qc_new.b = r;
            qc_new.lce_b = lce_r;
        } else {
            qc_new.b = l;
            qc_new.lce_b = lce_l;
        }
    } else {
        qc_new.b = b;
        qc_new.lce_b = lce_b;
    }

    if (lce_e < len) {
        sidx_t l = e_min;
        sidx_t r = e_max;
        pos_t lce_l = lce_e_min;
        pos_t lce_r = lce_e_max;

        sidx_t m;
        pos_t pos_m, lce_m;

        while (r - l > 1) {
            m = l + (r - l) / 2;
            pos_m = S[SXA<dir>(m)];

            lce_m = lce_offs<dir>(
                pos_pat, pos_m,
                std::min<pos_t>(lce_l, lce_r),
                len
            );

            if (lce_m >= len) {
                l = m;
                lce_l = lce_m;
            } else {
                r = m;
                lce_r = lce_m;
            }
        }

        if (lce_r >= len) {
            qc_new.e = r;
            qc_new.lce_e = lce_r;
        } else {
            qc_new.e = l;
            qc_new.lce_e = lce_l;
        }
    } else {
        qc_new.e = e;
        qc_new.lce_e = lce_e;
    }

    return true;
}

template <typename pos_t, typename sidx_t, typename lce_r_t>
template <direction dir>
sample_index<pos_t, sidx_t, lce_r_t>::query_context
sample_index<pos_t, sidx_t, lce_r_t>::interpolate(
    const query_context& qc_short,
    const query_context& qc_long,
    pos_t pos_pat, pos_t len
) const {
    if (qc_short.match_length() >= len) {
        return qc_short;
    }

    sidx_t l = qc_short.b;
    sidx_t r = qc_long.b;
    sidx_t m;

    pos_t lce_l = qc_short.lce_b;
    pos_t lce_r = qc_long.lce_b;
    pos_t lce_m;

    query_context qc_ret;
    pos_t pos_m;

    if constexpr (dir == LEFT) {
        lce_l = lce_offs<dir>(
            pos_pat, S[SXA<dir>(l)],
            lce_l,
            len
        );
    }

    while (r - l > 1) {
        m = l + (r - l) / 2;
        pos_m = S[SXA<dir>(m)];
        lce_m = lce_offs<dir>(
            pos_pat, pos_m,
            std::min<pos_t>(lce_l, lce_r),
            len
        );

        if (lce_m < len) {
            l = m;
            lce_l = lce_m;
        } else {
            r = m;
            lce_r = lce_m;
        }
    }

    if (lce_l < len) {
        qc_ret.b = r;
        qc_ret.lce_b = lce_r;
    } else {
        qc_ret.b = l;
        qc_ret.lce_b = lce_l;
    }

    l = qc_long.e;
    r = qc_short.e;

    lce_l = qc_long.lce_e;
    lce_r = qc_short.lce_e;

    if constexpr (dir == LEFT) {
        lce_r = lce_offs<dir>(
            pos_pat, S[SXA<dir>(r)],
            lce_r, len);
    }

    while (r - l > 1) {
        m = l + (r - l) / 2;
        pos_m = S[SXA<dir>(m)];

        lce_m = lce_offs<dir>(
            pos_pat, pos_m,
            std::min<pos_t>(lce_l, lce_r),
            len
        );

        if (lce_m < len) {
            r = m;
            lce_r = lce_m;
        } else {
            l = m;
            lce_l = lce_m;
        }
    }

    if (lce_r < len) {
        qc_ret.e = l;
        qc_ret.lce_e = lce_l;
    } else {
        qc_ret.e = r;
        qc_ret.lce_e = lce_r;
    }

    return qc_ret;
}