#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
template <typename output_fnc_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::
    transform_to_exact_without_samples(output_fnc_t output)
{
    if (log) {
        std::cout << "computing the exact factorization" << std::flush;
    }

    num_fact = 0;

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

        sidx_t x_r = 0;
        sidx_t x_c = bin_search_min_geq<pos_t, sidx_t>(
            b, 0, c - 1, [&](sidx_t x) { return C[x]; });
        pos_t num_phr_thr = 0;
        std::ofstream fact_ofile;
        if (p > 1) fact_ofile.open(fact_file_name + "_" + std::to_string(i_p));
        std::ostream_iterator<factor> fact_it(fact_ofile);

        for (pos_t i = b; i < e;) {
            while (beg_nxt_aprx_phr <= i) {
                f_aprx = *aprx_it++;
                beg_nxt_aprx_phr += f_aprx.length();
            }

            factor f = f_aprx;
            pos_t max_j = std::min<pos_t>(e, i + delta);

            if (f.len != 0) {
                pos_t cut_left = f.len - (beg_nxt_aprx_phr - i);
                f.len -= cut_left;
                f.src += cut_left;
            }

            if constexpr (range_ds_t<sidx_t>::is_dynamic()) {
                insert_points(x_r, i);
                find_close_sources(f, i, e);
            }

            #if defined(GEN_RANGE_QUERIES)
            insert_points(x_r, i);
            #endif

            for (pos_t j = i; j < max_j; j++) {
                pos_t lce_l = (j - i) + 1;
                query_context_t qc_left = idx_C.query();

                if (idx_C.extend_left(qc_left, j, lce_l, false)) {
                    query_context_t qc_right = idx_C.query();
                    query_context_t qc_right_nxt = idx_C.query();

                    pos_t lce_r_min = f.len < j - i ? 0 : (i + f.len - j);
                    pos_t lce_r_max = e - j;

                    exp_search_max_geq<bool, pos_t, RIGHT>(
                        true, lce_r_min, lce_r_max,
                        [&](pos_t lce_r_tmp) {
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

            #ifndef NDEBUG
            assert((f.len == 0 && (char) f.src == T[i]) || f.src < i);
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