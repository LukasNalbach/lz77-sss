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
    improve_factor_without_interval_samples(uint64_t i, uint64_t e, uint64_t& x_c, factor& f)
{
    uint64_t max_j = std::min<uint64_t>(e, i + delta);

    for (uint64_t j = i; j < max_j; j++) {
        if (f.len >= e - i) [[unlikely]] break;
        uint64_t lce_l = (j - i) + 1;
        query_ctx_t qc_left = idx_C.query();

        if (idx_C.extend_left(qc_left, j, lce_l, interval_samples::skip)) {
            query_ctx_t qc_right = idx_C.query();
            query_ctx_t qc_right_nxt = idx_C.query();

            uint64_t lce_r_min = f.len < j - i ? 0 : (i + f.len - j);
            uint64_t lce_r_max = e - j;

            exp_search_max_geq<bool, uint64_t, RIGHT>(
                true, lce_r_min, lce_r_max,
                [&](uint64_t lce_r_tmp) {
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
                        if (intersect(qc_left.interval(),
                                qc_right_tmp.interval(), i, j,
                                lce_l, lce_r_tmp, x_c, f)) {
                            qc_right = qc_right_tmp;
                            return true;
                        } else {
                            qc_right_nxt = qc_right_tmp;
                        }
                    }

                    return false;
                });
        }
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::exact_transformer::
    transform_to_exact_without_interval_samples(factor_sink& output)
{
    log_phase_begin(log, "computing the exact factorization");
    phase_progress progress(log, n);

    num_fact = 0;

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
    for (uint64_t sect = 0; sect < num_par_sect; sect++) {
        uint64_t b = par_sect[sect].beg;
        uint64_t e = par_sect[sect + 1].beg;

        direct_ifstream aprx_ifile(aprx_file_name);
        aprx_ifile.seekg(par_sect[sect].first_aprx * (src_bytes + len_bytes), std::ios::beg);
        factor f_aprx;
        f_aprx.read(aprx_ifile, src_bytes, len_bytes);
        uint64_t aprx_end = b + f_aprx.text_len();

        uint64_t x_r = 0;
        uint64_t x_c = bin_search_min_geq<uint64_t, uint64_t>(
            b, 0, c - 1, [&](uint64_t x) { return C[x]; });
        uint64_t num_fact_sect = 0;
        direct_ofstream fact_ofile;
        if (p > 1) fact_ofile.open(fact_file_name + "_" + std::to_string(sect));

        for (uint64_t i = b; i < e;) {
            while (aprx_end <= i) {
                f_aprx.read(aprx_ifile, src_bytes, len_bytes);
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

            if (write_queries) insert_points_before(x_r, i);

            improve_factor_without_interval_samples(i, e, x_c, f);

            #ifndef NDEBUG
            assert((f.is_literal() && f.src == T.to_char(T[i])) || f.src < i);
            assert(f.len <= e - i);

            for (uint64_t j = 0; j < f.len; j++) {
                assert(T[f.src + j] == T[i + j]);
            }
            #endif

            if (p == 1) output(f);
            else f.write(fact_ofile, src_bytes, len_bytes);
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
        combine_factorizations(output, [&](uint64_t i) {
            factor f = aprx_factor_at(i);
            uint64_t x_c = bin_search_min_geq<uint64_t, uint64_t>(
                i, 0, c - 1, [&](uint64_t x) { return C[x]; });
            improve_factor_without_interval_samples(i, n, x_c, f);
            return f;
        });
    }
}