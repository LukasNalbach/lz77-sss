#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::
extend_right_with_samples(
    const sxa_interval_t& spa_iv,
    pos_t i, pos_t j, sidx_t& x_c, factor& f
) {
    const rk61_substring& rks = idx_C.rabin_karp_substring();
    const std::vector<pos_t>& smpl_lens_right = idx_C.sampled_pattern_lengths_right();
    pos_t num_smpl_lens_right = idx_C.num_sampled_pattern_lengths_right();
    pos_t lce_r_min = f.len < j - i ? 0 : (i + f.len - j);
    pos_t lce_l = (j - i) + 1;

    int16_t x_min = bin_search_max_leq<pos_t, int16_t>(
        lce_r_min, 0, num_smpl_lens_right - 1, [&](int16_t x){
            return smpl_lens_right[x];
    });

    int16_t x_max = bin_search_max_leq<pos_t, int16_t>(
        n - j, x_min, num_smpl_lens_right - 1, [&](int16_t x){
            return smpl_lens_right[x];
    });

    sxa_interval_t ssa_iv;
    sxa_interval_t ssa_iv_nxt;
    pos_t lce_r = 0;
    pos_t lce_r_nxt = 0;
    uint64_t fp_right = 0;
    
    int16_t x_res = exp_search_max_geq<bool, int16_t, RIGHT>(true, x_min - 1, x_max, 1, [&](int16_t x){
        time_point_t t2;
        if (log) t2 = now();
        pos_t lce_r_tmp = smpl_lens_right[x];
        pos_t len_add = lce_r_tmp - lce_r;
        uint64_t fp_add = rks.substring(j + lce_r, len_add);
        uint64_t fp_tmp = rks.concat(fp_right, fp_add, len_add);
        auto [ssa_iv_tmp, result] = idx_C.ssa_interval(x, j, fp_tmp);

        if (result) {
            if (log) time_extend_right += time_diff_ns(t2);

            if (intersect(spa_iv, ssa_iv_tmp, i, j, lce_l, lce_r_tmp, x_c, f)) {
                ssa_iv = ssa_iv_tmp;
                fp_right = fp_tmp;
                lce_r = lce_r_tmp;
                return true;
            } else {
                ssa_iv_nxt = ssa_iv_tmp;
                lce_r_nxt = lce_r_tmp;
            }
        } else if (log) {
            time_extend_right += time_diff_ns(t2);
        }

        return false;
    });

    if (x_res < x_min || (x_res < x_max && smpl_lens_right[x_res + 1] < lce_r_min)) {
        return;
    }

    query_context_t qc_right = idx_C.query_right(ssa_iv, j, lce_r);
    query_context_t qc_right_nxt = lce_r_nxt == 0 ? idx_C.query() :
        idx_C.query_right(ssa_iv_nxt, j, lce_r_nxt);
    assert(x_res < 0 || qc_right.match_length() >= smpl_lens_right[x_res]);
    pos_t lce_r_max = lce_r_nxt == 0 ? n - j : (lce_r_nxt  - 1);

    exp_search_max_geq<bool, pos_t, RIGHT>(true, lce_r, lce_r_max, 1, [&](pos_t lce_r_tmp){
        query_context_t qc_right_tmp;
        time_point_t t2;
        bool result = true;
        if (log) t2 = now();
        
        if (qc_right_nxt.match_length() == 0) {
            result = idx_C.extend_right(
                qc_right, qc_right_tmp, j, lce_r_tmp, false);
        } else {
            qc_right_tmp = idx_C.interpolate_right(
                qc_right, qc_right_nxt, j, lce_r_tmp);
        }

        if (result) {
            if (log) time_extend_right += time_diff_ns(t2);

            if (intersect(spa_iv, qc_right_tmp.interval(),
                i, j, lce_l, lce_r_tmp, x_c, f)
            ) {
                qc_right = qc_right_tmp;
                return true;
            } else {
                qc_right_nxt = qc_right_tmp;
            }
        } else if (log) {
            time_extend_right += time_diff_ns(t2);
        }

        return false;
    });
}

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::
transform_to_exact_with_samples(out_it_t& out_it) {
    if (log) {
        std::cout << "computing the exact factorization" << std::flush;
    }
    
    std::vector<uint64_t> fp_left(delta);
    const rk61_substring& rks = idx_C.rabin_karp_substring();
    const std::vector<pos_t>& smpl_lens_left = idx_C.sampled_pattern_lengths_left();
    std::vector<uint8_t> is_smpld_left(delta + 1, 0);

    for (pos_t len : smpl_lens_left) {
        is_smpld_left[len] = 1;
    }

    sidx_t x_c = 0;
    num_phr = 0;

    for (pos_t i = 0; i < n;) {
        factor f {.src = char_to_uchar(T[i]), .len = 0};

        if constexpr (is_dynamic<range_ds_t>()) {
            insert_points(x_c, i);
            handle_close_sources(f, i);
        }

        pos_t max_k = std::min<pos_t>(delta, n - i);
        fp_left[0] = char_to_uchar(T[i]);

        for (pos_t k = 1; k < max_k; k++) {
            fp_left[k] = rks.push(fp_left[k - 1], T[i + k]);
        }

        for (pos_t x = 0; x < smpl_lens_left.size(); x++) {
            pos_t k = smpl_lens_left[x] - 1;
            pos_t j = i + k;
            time_point_t t1;
            if (log) t1 = now();
            auto [spa_iv, result] = idx_C.spa_interval(x, j, fp_left[k]);

            if (result) {
                if (log) time_extend_left += time_diff_ns(t1);
                extend_right_with_samples(spa_iv, i, j, x_c, f);
            } else if (log) {
                time_extend_left += time_diff_ns(t1);
            }
        }

        for (pos_t k = 2; k < max_k; k++) {
            pos_t lce_l = k + 1;
            if (is_smpld_left[lce_l]) continue;
            pos_t j = i + k;
            time_point_t t1;
            if (log) t1 = now();
            query_context_t qc_left = idx_C.query();

            if (idx_C.extend_left(qc_left, j, lce_l)) {
                if (log) time_extend_left += time_diff_ns(t1);
                extend_right_with_samples(qc_left.interval(), i, j, x_c, f);
            } else if (log) {
                time_extend_left += time_diff_ns(t1);
            }
        }

        #ifndef NDEBUG
        assert((f.len == 0 && (char) f.src == T[i]) || f.src < i);

        for (pos_t j = 0; j < f.len; j++) {
            assert(T[f.src + j] == T[i + j]);
        }
        #endif

        *out_it++ = f;
        i += std::max<pos_t>(1, f.len);
        num_phr++;
    }

    if (log) {
        time = log_runtime(time);
    }
}