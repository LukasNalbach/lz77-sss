#pragma once

#include <lz77_sss/data_structures/ring_buffer.hpp>

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
bool lz77_sss<pos_t>::factorizer<tau, out_it_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::
intersect(
    const sxa_interval_t& spa_iv, const sxa_interval_t& ssa_iv,
    pos_t i, pos_t j, pos_t lce_l, pos_t lce_r, sidx_t& x_c, factor& f
) {
    time_point_t t3;
    point_t p;
    bool result = false;
    if (log) t3 = now();
    adjust_sample_index(x_c, j);

    pos_t spa_rng = spa_iv.e - spa_iv.b + 1;
    pos_t ssa_rng = ssa_iv.e - ssa_iv.b + 1;
    
    if (spa_rng <= range_scan_threshold && spa_rng <= ssa_rng) {
        for (sidx_t x = spa_iv.b; x <= spa_iv.e; x++) {
            if (idx_C.spa(x) < x_c &&
                ssa_iv.b <= PS[x] && PS[x] <= ssa_iv.e
            ) {
                p.x = x;
                p.y = PS[x];
                result = true;
                break;
            }
        }
    } else if (ssa_rng <= range_scan_threshold) {
        for (sidx_t y = ssa_iv.b; y <= ssa_iv.e; y++) {
            if (idx_C.ssa(y) < x_c &&
                spa_iv.b <= SP[y] && SP[y] <= spa_iv.e
            ) {
                p.x = SP[y];
                p.y = y;
                result = true;
                break;
            }
        }
    } else if constexpr (is_static<range_ds_t>()) {
        std::tie(p, result) = R.lighter_point_in_range(
            x_c, spa_iv.b, spa_iv.e, ssa_iv.b, ssa_iv.e
        );
    } else {
        std::tie(p, result) = R.point_in_range(
            spa_iv.b, spa_iv.e, ssa_iv.b, ssa_iv.e
        );
    }

    if (result) {
        #ifndef NDEBUG
        assert(spa_iv.b      <= p.x && p.x <= spa_iv.e);
        assert(ssa_iv.b <= p.y && p.y <= ssa_iv.e);
        #endif

        pos_t lce = lce_l + lce_r - 1;

        if (lce > f.len) {
            f.len = lce;
            f.src = C[idx_C.ssa(p.y)] - lce_l + 1;

            #ifndef NDEBUG
            assert(f.src < i);

            for (pos_t x = 0; x < f.len; x++) {
                assert(T[i + x] == T[f.src + x]);
            }
            #endif
        }
    }

    if (log) {
        num_range_queries++;
        spa_range_sum += spa_iv.e - spa_iv.b + 1;
        ssa_range_sum += ssa_iv.e - ssa_iv.b + 1;
        time_range_queries += time_diff_ns(t3);
    }
        
    return result;
}

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, out_it_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::
extend_right(
    const sample_index<pos_t, sidx_t, lce_t>::sxa_interval_t& spa_iv, pos_t i, pos_t j, sidx_t& x_c, factor& f
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
    pos_t last_lce_r = 0;
    uint64_t fp_right = 0;
    
    int16_t x_res = exp_search_max_geq<bool, int16_t, RIGHT>(true, x_min - 1, x_max, 1, [&](int16_t x){
        time_point_t t2;
        if (log) t2 = now();
        pos_t lce_r = smpl_lens_right[x];
        pos_t len_add = lce_r - last_lce_r;
        uint64_t fp_add = rks.substring(j + last_lce_r, len_add);
        uint64_t fp_right_tmp = rks.concat(fp_right, fp_add, len_add);
        auto [ssa_iv_tmp, result] = idx_C.ssa_interval(x, j, fp_right_tmp);

        if (result) {
            if (log) time_extend_right += time_diff_ns(t2);

            if (intersect(spa_iv, ssa_iv_tmp, i, j, lce_l, lce_r, x_c, f)) {
                ssa_iv = ssa_iv_tmp;
                fp_right = fp_right_tmp;
                last_lce_r = lce_r;
                return true;
            }
        } else if (log) {
            time_extend_right += time_diff_ns(t2);
        }

        return false;
    });

    if (x_res < x_min || (x_res < x_max && smpl_lens_right[x_res + 1] < lce_r_min)) {
        return;
    }

    query_context_t qc_right = idx_C.query_right(ssa_iv, j, last_lce_r);

    exp_search_max_geq<bool, pos_t, RIGHT>(true, last_lce_r, n - j, 1, [&](pos_t lce_r){
        query_context_t qc_right_tmp;
        time_point_t t2;
        if (log) t2 = now();

        if (idx_C.extend_right(qc_right, qc_right_tmp, j, lce_r)) {
            if (log) time_extend_right += time_diff_ns(t2);

            if (intersect(spa_iv, qc_right_tmp.interval(), i, j, lce_l, lce_r, x_c, f)) {
                qc_right = qc_right_tmp;
                return true;
            }
        } else if (log) {
            time_extend_right += time_diff_ns(t2);
        }

        return false;
    });
}

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, out_it_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::
transform_to_exact_optimized(out_it_t& out_it) {
    if (log) {
        std::cout << "computing the exact factorization" << std::flush;
    }
    
    std::vector<uint64_t> fp_left(delta);
    const rk61_substring& rks = idx_C.rabin_karp_substring();
    const std::vector<pos_t>& smpl_lens_left = idx_C.sampled_pattern_lengths_left();
    std::vector<uint8_t> is_smpld_left(delta, 0);

    for (pos_t len : smpl_lens_left) {
        is_smpld_left[len] = 1;
    }

    sidx_t x_c = 0;
    num_phr = 0;

    for (pos_t i = 0; i < n;) {
        factor f {.src = char_to_uchar(T[i]), .len = 0};

        if constexpr (is_dynamic<range_ds_t>()) {
            while (x_c < c && C[x_c] < i) {
                time_point_t t0;
                if (log) t0 = now();
                R.insert(P[x_c]);
                if (log) time_insert_points += time_diff_ns(t0);
                x_c++;
            }

            pos_t min_j = i <= delta ? 0 : (i - delta);
            time_point_t t5;
            if (log) t5 = now();

            for (pos_t j = min_j; j < i; j++) {
                pos_t lce = LCE_R(j, i);
                
                if (lce > f.len) {
                    f.src = j;
                    f.len = lce;
                }
            }

            if (log) time_close_sources += time_diff_ns(t5);
        }

        pos_t max_k = std::min<pos_t>(delta, n - i);
        fp_left[0] = T[i];

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
                extend_right(spa_iv, i, j, x_c, f);
            } else if (log) {
                time_extend_left += time_diff_ns(t1);
            }
        }

        for (pos_t k = 0; k < max_k; k++) {
            pos_t lce_l = k + 1;
            if (is_smpld_left[lce_l]) continue;
            pos_t j = i + k;
            time_point_t t1;
            if (log) t1 = now();
            query_context_t qc_left = idx_C.query();

            if (idx_C.extend_left(qc_left, j, lce_l)) {
                if (log) time_extend_left += time_diff_ns(t1);
                extend_right(qc_left.interval(), i, j, x_c, f);
            } else if (log) {
                time_extend_left += time_diff_ns(t1);
            }
        }

        *out_it++ = f;
        i += std::max<pos_t>(1, f.len);
        num_phr++;
    }

    if (log) {
        time = log_runtime(time);
    }
}