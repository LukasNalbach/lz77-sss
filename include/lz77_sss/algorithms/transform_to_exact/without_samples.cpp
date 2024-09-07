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

    sidx_t x_c = 0;
    num_phr = 0;

    for (pos_t i = 0; i < n;) {
        factor f {.src = char_to_uchar(T[i]), .len = 0};

        if constexpr (range_ds_t<sidx_t>::is_dynamic()) {
            insert_points(x_c, i);
            handle_close_sources(f, i);
        }

        pos_t max_j = std::min<pos_t>(n, i + delta);

        for (pos_t j = i; j < max_j; j++) {
            pos_t lce_l = (j - i) + 1;

            query_context_t qc_left = idx_C.query();
            time_point_t t1;
            if (log) t1 = now();

            if (idx_C.extend_left(qc_left, j, lce_l, false)) {
                if (log) time_extend_left += time_diff_ns(t1);
                query_context_t qc_right = idx_C.query();
                query_context_t qc_right_nxt = idx_C.query();
                pos_t lce_r_min = f.len < j - i ? 0 : (i + f.len - j);
                pos_t lce_r_max = n - j;

                exp_search_max_geq<bool, pos_t, RIGHT>(
                    true, lce_r_min, lce_r_max,
                    1, [&](pos_t lce_r_tmp)
                {
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

                        if (intersect(qc_left.interval(),
                            qc_right_tmp.interval(), i, j,
                            lce_l, lce_r_tmp, x_c, f)
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