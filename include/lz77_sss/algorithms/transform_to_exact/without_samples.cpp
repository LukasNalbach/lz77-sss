#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::
transform_to_exact_without_samples(out_it_t& out_it) {
    if (log) {
        std::cout << "computing the exact factorization" << std::flush;
    }

    num_phr = 0;

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b = start_thr[i_p];
        pos_t e = start_thr[i_p + 1];

        sidx_t x_c = bin_search_min_geq<pos_t, sidx_t>(
            b, 0, c - 1, [&](sidx_t x) {return C[x];});
        sidx_t x_r = 0;
        pos_t num_phr_thr = 0;
        std::ofstream fact_ofile;
        if (p > 1) fact_ofile.open(
            fact_file_name + "_" + std::to_string(i_p));
        std::ostream_iterator<factor> fact_it(fact_ofile);

        for (pos_t i = b; i < e;) {
            factor f {.src = char_to_uchar(T[i]), .len = 0};
            pos_t max_j = std::min<pos_t>(e, i + delta);

            if constexpr (range_ds_t<sidx_t>::is_dynamic()) {
                insert_points(x_r, i);
                find_close_sources(f, i, e);
            }

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
                        [&](pos_t lce_r_tmp)
                    {
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
                                lce_l, lce_r_tmp, x_c, f)
                            ) {
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

            if (p == 1) *out_it++ = f;
            else *fact_it++ = f;
            i += std::max<pos_t>(1, f.len);
            num_phr_thr++;
        }

        #pragma omp critical
        {num_phr += num_phr_thr;}
    }

    if (p > 1) {
        combine_factorizations(out_it);
    }

    if (log) {
        time = log_runtime(time);
    }
}