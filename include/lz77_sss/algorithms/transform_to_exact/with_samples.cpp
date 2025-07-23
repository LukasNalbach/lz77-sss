#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::
    extend_right_with_samples(
        const interval_t& pa_c_iv,
        pos_t i, pos_t j, pos_t e, sidx_t& x_c, factor& f)
{
    const auto& rks = idx_C.rab_karp_substr();
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

    interval_t sa_c_iv;
    interval_t sa_c_iv_nxt { .b = 1, .e = 0 };
    pos_t lce_r = 0;
    pos_t lce_r_nxt = 0;
    uint64_t fp_right = 0;

    int16_t x_res = exp_search_max_geq<bool, int16_t, RIGHT>(true, x_min - 1, x_max, [&](int16_t x) {
        pos_t lce_r_tmp = smpl_lens_right[x];
        pos_t len_add = lce_r_tmp - lce_r;
        uint64_t fp_add = rks.substring(j + lce_r, len_add);
        uint64_t fp_tmp = rks.concat(fp_right, fp_add, len_add);
        auto [sa_c_iv_tmp, result] = idx_C.sa_s_interval(x, j, fp_tmp);

        if (result) {
            if (intersect(pa_c_iv, sa_c_iv_tmp, i, j, lce_l, lce_r_tmp, x_c, f)) {
                sa_c_iv = sa_c_iv_tmp;
                fp_right = fp_tmp;
                lce_r = lce_r_tmp;
                return true;
            } else {
                sa_c_iv_nxt = sa_c_iv_tmp;
            }
        }

        lce_r_nxt = lce_r_tmp;
        return false;
    });

    if (x_res < x_min || (x_res < x_max && smpl_lens_right[x_res + 1] < lce_r_min)) {
        return;
    }

    query_context_t qc_right = idx_C.query_right(sa_c_iv, j, lce_r);
    query_context_t qc_right_nxt = sa_c_iv_nxt.empty() ? idx_C.query() : idx_C.query_right(sa_c_iv_nxt, j, lce_r_nxt);
    assert(x_res < 0 || qc_right.match_length() >= smpl_lens_right[x_res]);
    pos_t lce_r_max = lce_r_nxt == 0 ? e - j : (lce_r_nxt - 1);

    auto fnc = [&](auto lce_r_tmp) {
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
            if (intersect(pa_c_iv, qc_right_tmp.interval(),
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
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
template <typename output_fnc_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::
    transform_to_exact_with_samples(output_fnc_t output)
{
    if (log) {
        std::cout << "computing the exact factorization" << std::flush;
    }

    const auto& rks = idx_C.rab_karp_substr();
    const std::vector<pos_t>& smpl_lens_left = idx_C.sampled_pattern_lengths_left();
    std::vector<uint8_t> is_smpld_left(delta + 1, 0);
    num_fact = 0;

    for (pos_t len : smpl_lens_left) {
        if (len <= delta) {
            is_smpld_left[len] = 1;
        }
    }

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b = par_sect[i_p].beg;
        pos_t e = par_sect[i_p + 1].beg;

        std::ifstream aprx_ifile(aprx_file_name);
        aprx_ifile.seekg(par_sect[i_p].phr_idx *
            lz77_sss<pos_t>::factor::size_of(), std::ios::beg);
        std::istream_iterator<factor> aprx_it(aprx_ifile);
        lz77_sss<pos_t>::factor f_aprx = *aprx_it++;
        pos_t beg_nxt_aprx_phr = b + f_aprx.length();

        sidx_t x_c = bin_search_min_geq<pos_t, sidx_t>(
            b, 0, c - 1, [&](sidx_t x) { return C[x]; });
        sidx_t x_r = 0;
        pos_t num_phr_thr = 0;
        std::ofstream fact_ofile;
        if (p > 1) fact_ofile.open(fact_file_name + "_" + std::to_string(i_p));
        std::ostream_iterator<factor> fact_it(fact_ofile);
        std::vector<uint64_t> fp_left(delta);

        for (pos_t i = b; i < e;) {
            while (beg_nxt_aprx_phr <= i) {
                f_aprx = *aprx_it++;
                beg_nxt_aprx_phr += f_aprx.length();
            }

            factor f = f_aprx;

            if (f.len != 0) {
                pos_t cut_left = f.len - (beg_nxt_aprx_phr - i);
                f.len -= cut_left;
                f.src += cut_left;
            }

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
                auto [pa_c_iv, result] = idx_C.pa_s_interval(x, j, fp_left[k]);
                if (result) extend_right_with_samples(pa_c_iv, i, j, e, x_c, f);
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
            assert((f.len == 0 && f.src == char_to_uchar(T[i])) || f.src < i);
            assert(f.len <= e - i);

            for (pos_t j = 0; j < f.len; j++) {
                assert(T[f.src + j] == T[i + j]);
            }
            #endif

            if (p == 1) output(f);
            else *fact_it++ = f;
            i += f.length();
            num_phr_thr++;
        }

        #pragma omp critical
        {
            num_fact += num_phr_thr;
        }
    }

    if (p > 1) {
        combine_factorizations(output);
    }
}