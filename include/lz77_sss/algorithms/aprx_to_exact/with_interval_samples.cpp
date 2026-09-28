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

#include <lz77_sss/lz77_sss.hpp>

template <typename text_t>
void lz77_sss::factorizer<text_t>::exact_transformer::
    extend_right_with_interval_samples(
        const interval_t& pa_c_iv,
        uint64_t i, uint64_t j, uint64_t e, uint64_t& x_c, factor& f)
{
    const auto& rks = idx_C.rks();
    const std::vector<uint64_t>& smpl_lens_right = idx_C.sampled_pattern_lengths_right();
    uint64_t num_smpl_lens_right = idx_C.num_sampled_pattern_lengths_right();
    uint64_t lce_r_min = f.len < j - i ? 0 : (i + f.len - j);
    uint64_t lce_l = (j - i) + 1;

    int16_t x_min = bin_search_max_leq<uint64_t, int16_t>(
        lce_r_min, 0, num_smpl_lens_right - 1, [&](int16_t x) {
            return smpl_lens_right[x];
        });

    int16_t x_max = bin_search_max_leq<uint64_t, int16_t>(
        e - j, x_min, num_smpl_lens_right - 1, [&](int16_t x) {
            return smpl_lens_right[x];
        });

    interval_t sa_c_iv;
    interval_t sa_c_iv_nxt { .b = 1, .e = 0 };
    uint64_t lce_r = 0;
    uint64_t lce_r_nxt = 0;
    uint64_t fp_right = 0;

    int16_t x_res = exp_search_max_geq<bool, int16_t, RIGHT>(true, x_min - 1, x_max, [&](int16_t x) {
        uint64_t lce_r_tmp = smpl_lens_right[x];
        uint64_t len_add = lce_r_tmp - lce_r;
        uint64_t fp_add = rks.substr_fp(j + lce_r, len_add);
        uint64_t fp_tmp = rks.concat(fp_right, fp_add, len_add);
        auto [sa_c_iv_tmp, result] = idx_C.sa_interval(x, j, fp_tmp);

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

    query_ctx_t qc_right = idx_C.query_right(sa_c_iv, j, lce_r);
    query_ctx_t qc_right_nxt = sa_c_iv_nxt.empty() ? idx_C.query() : idx_C.query_right(sa_c_iv_nxt, j, lce_r_nxt);
    assert(x_res < 0 || qc_right.match_length() >= smpl_lens_right[x_res]);
    uint64_t lce_r_max = lce_r_nxt == 0 ? e - j : (lce_r_nxt - 1);

    auto fnc = [&](auto lce_r_tmp) {
        query_ctx_t qc_right_tmp;
        bool result = true;

        if (qc_right_nxt.match_length() == 0) {
            result = idx_C.extend_right(
                qc_right, qc_right_tmp, j, lce_r_tmp, interval_samples::skip);
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
        exp_search_max_geq<bool, uint64_t, RIGHT>(true, lce_r, lce_r_max, fnc);
    } else {
        bin_search_max_geq<bool, uint64_t>(true, lce_r, lce_r_max, fnc);
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::exact_transformer::
    improve_factor_with_interval_samples(uint64_t i, uint64_t e, uint64_t& x_c, std::vector<uint64_t>& fp_left, factor& f)
{
    const auto& rks = idx_C.rks();
    const std::vector<uint64_t>& smpl_lens_left = idx_C.sampled_pattern_lengths_left();
    uint64_t max_k = std::min<uint64_t>(delta, e - i);
    fp_left[0] = T[i];

    for (uint64_t k = 1; k < max_k; k++) {
        fp_left[k] = rks.push(fp_left[k - 1], T[i + k]);
    }

    for (uint64_t x = 0; x < smpl_lens_left.size(); x++) {
        uint64_t k = smpl_lens_left[x] - 1;
        if (k >= max_k) break;
        uint64_t j = i + k;
        auto [pa_c_iv, result] = idx_C.pa_interval(x, j, fp_left[k]);
        if (result) extend_right_with_interval_samples(pa_c_iv, i, j, e, x_c, f);
    }

    for (uint64_t k = 2; k < max_k; k++) {
        if (f.len >= e - i) [[unlikely]] break;
        uint64_t lce_l = k + 1;
        if (is_smpl_len_left[lce_l]) continue;
        uint64_t j = i + k;
        query_ctx_t qc_left = idx_C.query();

        if (idx_C.extend_left(qc_left, j, lce_l)) {
            extend_right_with_interval_samples(qc_left.interval(), i, j, e, x_c, f);
        }
    }

    if (f.len > e - i) [[unlikely]] {
        f.len = e - i;
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::exact_transformer::
    transform_to_exact_with_interval_samples(factor_sink& output)
{
    log_phase_begin(log, "computing the exact factorization");
    phase_progress progress(log, n);

    const std::vector<uint64_t>& smpl_lens_left = idx_C.sampled_pattern_lengths_left();
    is_smpl_len_left.assign(delta + 1, 0);
    num_fact = 0;

    for (uint64_t len : smpl_lens_left) {
        if (len <= delta) {
            is_smpl_len_left[len] = 1;
        }
    }

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t sect = 0; sect < num_par_sect; sect++) {
        uint64_t b = par_sect[sect].beg;
        uint64_t e = par_sect[sect + 1].beg;

        direct_ifstream aprx_ifile(aprx_file_name);
        aprx_ifile.seekg(par_sect[sect].first_aprx *
            lz77_sss::factor::size_of(), std::ios::beg);
        std::istream_iterator<factor> aprx_it(aprx_ifile);
        lz77_sss::factor f_aprx = *aprx_it++;
        uint64_t aprx_end = b + f_aprx.text_len();

        uint64_t x_c = bin_search_min_geq<uint64_t, uint64_t>(
            b, 0, c - 1, [&](uint64_t x) { return C[x]; });
        uint64_t x_r = 0;
        uint64_t num_fact_sect = 0;
        direct_ofstream fact_ofile;
        if (p > 1) fact_ofile.open(fact_file_name + "_" + std::to_string(sect));
        std::ostream_iterator<factor> fact_it(fact_ofile);
        std::vector<uint64_t> fp_left(delta);

        for (uint64_t i = b; i < e;) {
            while (aprx_end <= i) {
                f_aprx = *aprx_it++;
                aprx_end += f_aprx.text_len();
            }

            factor f = f_aprx;

            if (!f.is_literal()) {
                uint64_t cut_left = f.len - (aprx_end - i);
                f.len = f.len - cut_left;
                f.src = f.src + cut_left;
            }

            if (kind.is_dynamic()) {
                insert_points_before(x_r, i);
                find_close_sources(f, i, e);
            }

            improve_factor_with_interval_samples(i, e, x_c, fp_left, f);

            #ifndef NDEBUG
            assert((f.is_literal() && f.src == T.to_char(T[i])) || f.src < i);
            assert(f.len <= e - i);

            for (uint64_t j = 0; j < f.len; j++) {
                assert(T[f.src + j] == T[i + j]);
            }
            #endif

            if (p == 1) output(f);
            else *fact_it++ = f;
            i += f.text_len();
            progress.advance(f.text_len());
            num_fact_sect++;
        }

        #pragma omp critical
        {
            num_fact += num_fact_sect;
        }
    }

    if (p > 1) {
        std::vector<uint64_t> fp_left(delta);

        combine_factorizations(output, [&](uint64_t i) {
            factor f = aprx_factor_at(i);
            uint64_t x_c = bin_search_min_geq<uint64_t, uint64_t>(
                i, 0, c - 1, [&](uint64_t x) { return C[x]; });
            improve_factor_with_interval_samples(i, n, x_c, fp_left, f);
            return f;
        });
    }
}