#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::
    extend_right_with_samples(
        const sxa_interval_t& spa_iv,
        pos_t i, pos_t j, pos_t e, sidx_t& x_c, factor& f)
{
    const rk61_substring& rks = idx_C.rabin_karp_substring();
    const std::vector<pos_t>& smpl_lens_right = idx_C.sampled_pattern_lengths_right();
    pos_t num_smpl_lens_right = idx_C.num_sampled_pattern_lengths_right();
    pos_t lce_r_min = f.len < j - i ? 0 : (i + f.len - j);
    pos_t lce_l = (j - i) + 1;

    int16_t x_min = bin_search_max_leq<pos_t, int16_t>(
        lce_r_min, 0, num_smpl_lens_right - 1, [&](int16_t x) {
            return smpl_lens_right[x];
        });

    int16_t x_max = bin_search_max_leq<pos_t, int16_t>(
        e - j, x_min, num_smpl_lens_right - 1, [&](int16_t x) {
            return smpl_lens_right[x];
        });

    sxa_interval_t ssa_iv;
    sxa_interval_t ssa_iv_nxt { .b = 1, .e = 0 };
    pos_t lce_r = 0;
    pos_t lce_r_nxt = 0;
    uint64_t fp_right = 0;

    int16_t x_res = exp_search_max_geq<bool, int16_t, RIGHT>(true, x_min - 1, x_max, [&](int16_t x) {
        pos_t lce_r_tmp = smpl_lens_right[x];
        pos_t len_add = lce_r_tmp - lce_r;
        uint64_t fp_add = rks.substring(j + lce_r, len_add);
        uint64_t fp_tmp = rks.concat(fp_right, fp_add, len_add);
        auto [ssa_iv_tmp, result] = idx_C.ssa_interval(x, j, fp_tmp);

        if (result) {
            if (intersect(spa_iv, ssa_iv_tmp, i, j, lce_l, lce_r_tmp, x_c, f)) {
                ssa_iv = ssa_iv_tmp;
                fp_right = fp_tmp;
                lce_r = lce_r_tmp;
                return true;
            } else {
                ssa_iv_nxt = ssa_iv_tmp;
            }
        }

        lce_r_nxt = lce_r_tmp;
        return false;
    });

    if (x_res < x_min || (x_res < x_max && smpl_lens_right[x_res + 1] < lce_r_min)) {
        return;
    }

    query_context_t qc_right = idx_C.query_right(ssa_iv, j, lce_r);
    query_context_t qc_right_nxt = ssa_iv_nxt.empty() ? idx_C.query() : idx_C.query_right(ssa_iv_nxt, j, lce_r_nxt);
    assert(x_res < 0 || qc_right.match_length() >= smpl_lens_right[x_res]);
    pos_t lce_r_max = lce_r_nxt == 0 ? e - j : (lce_r_nxt - 1);

    auto fnc = [&](pos_t lce_r_tmp) {
        query_context_t qc_right_tmp;
        bool result = true;

        if (qc_right_nxt.match_length() == 0) {
            result = idx_C.extend_right(
                qc_right, qc_right_tmp, j, lce_r_tmp, false);
        } else {
            qc_right_tmp = idx_C.interpolate_right(
                qc_right, qc_right_nxt, j, lce_r_tmp);
        }

        if (result) {
            if (intersect(spa_iv, qc_right_tmp.interval(),
                    i, j, lce_l, lce_r_tmp, x_c, f)) {
                qc_right = qc_right_tmp;
                return true;
            } else {
                qc_right_nxt = qc_right_tmp;
            }
        }

        return false;
    };

    if (lce_r_nxt == 0) {
        exp_search_max_geq<bool, pos_t, RIGHT>(true, lce_r, lce_r_max, fnc);
    } else {
        bin_search_max_geq<bool, pos_t>(true, lce_r, lce_r_max, fnc);
    }
}

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::
    transform_to_exact_with_samples(output_it_t& output)
{
    if (log) {
        std::cout << "computing the exact factorization" << std::flush;
    }

    const rk61_substring& rks = idx_C.rabin_karp_substring();
    const std::vector<pos_t>& smpl_lens_left = idx_C.sampled_pattern_lengths_left();
    std::vector<uint8_t> is_smpld_left(delta + 1, 0);
    num_phr = 0;

    for (pos_t len : smpl_lens_left) {
        is_smpld_left[len] = 1;
    }

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b = start_thr[i_p];
        pos_t e = start_thr[i_p + 1];

        sidx_t x_c = bin_search_min_geq<pos_t, sidx_t>(
            b, 0, c - 1, [&](sidx_t x) { return C[x]; });
        sidx_t x_r = 0;
        pos_t num_phr_thr = 0;
        std::ofstream fact_ofile;
        if (p > 1) fact_ofile.open(fact_file_name + "_" + std::to_string(i_p));
        std::ostream_iterator<factor> fact_it(fact_ofile);

        std::vector<uint64_t> fp_left(delta);

        for (pos_t i = b; i < e;) {
            factor f { .src = char_to_uchar(T[i]), .len = 0 };
            pos_t max_k = std::min<pos_t>(delta, e - i);
            fp_left[0] = char_to_uchar(T[i]);

            if constexpr (range_ds_t<sidx_t>::is_dynamic()) {
                insert_points(x_r, i);
                find_close_sources(f, i, e);
            }

            for (pos_t k = 1; k < max_k; k++) {
                fp_left[k] = rks.push(fp_left[k - 1], T[i + k]);
            }

            for (pos_t x = 0; x < smpl_lens_left.size(); x++) {
                pos_t k = smpl_lens_left[x] - 1;
                if (k >= max_k) break;
                pos_t j = i + k;
                auto [spa_iv, result] = idx_C.spa_interval(x, j, fp_left[k]);
                if (result) extend_right_with_samples(spa_iv, i, j, e, x_c, f);
            }

            for (pos_t k = 2; k < max_k; k++) {
                pos_t lce_l = k + 1;
                if (is_smpld_left[lce_l]) continue;
                pos_t j = i + k;
                query_context_t qc_left = idx_C.query();

                if (idx_C.extend_left(qc_left, j, lce_l)) {
                    extend_right_with_samples(qc_left.interval(), i, j, e, x_c, f);
                }
            }

            if (f.len > e - i) [[unlikely]] {
                f.len = e - i;
            }

            #ifndef NDEBUG
            assert((f.len == 0 && (char) f.src == T[i]) || f.src < i);
            assert(f.len <= e - i);

            for (pos_t j = 0; j < f.len; j++) {
                assert(T[f.src + j] == T[i + j]);
            }
            #endif

            if (p == 1) output(f);
            else *fact_it++ = f;
            i += std::max<pos_t>(1, f.len);
            num_phr_thr++;
        }

        #pragma omp critical
        {
            num_phr += num_phr_thr;
        }
    }

    if (p > 1) {
        combine_factorizations(output);
    }
}