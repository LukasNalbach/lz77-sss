#pragma once

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::transform_to_exact_naive(out_it_t& out_it) {
    if (log) {
        std::cout << "computing the exact factorization" << std::flush;
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
            
            pos_t min_j = i <= delta ? 0 : i - delta;
            time_point_t t5;
            if (log) t5 = now();

            for (pos_t j = min_j; j < i; j++) {
                pos_t lce_ji = LCE_R(j, i);
                
                if (lce_ji > f.len) {
                    f.src = j;
                    f.len = lce_ji;
                }
            }

            if (log) time_close_sources += time_diff_ns(t5);
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

                exp_search_max_geq<bool, pos_t, RIGHT>(true, 0, n - j, 1, [&](pos_t lce_r){
                    query_context_t qc_right_tmp;
                    time_point_t t2;
                    if (log) t2 = now();

                    if (idx_C.extend_right(qc_right, qc_right_tmp, j, lce_r, false)) {
                        time_point_t t3;
                        point_t sample_point;
                        bool result;

                        if (log) {
                            time_extend_right += time_diff_ns(t2);
                            t3 = now();
                        }

                        if constexpr (is_static<range_ds_t>()) {
                            adjust_sample_index(x_c, j);
                            std::tie(sample_point, result) = R.lighter_point_in_range(
                                x_c, qc_left.b, qc_left.e, qc_right_tmp.b, qc_right_tmp.e
                            );
                        } else {
                            std::tie(sample_point, result) = R.point_in_range(
                                qc_left.b, qc_left.e, qc_right_tmp.b, qc_right_tmp.e
                            );
                        }

                        if (result) {
                            assert(qc_left.b      <= sample_point.x && sample_point.x <= qc_left.e);
                            assert(qc_right_tmp.b <= sample_point.y && sample_point.y <= qc_right_tmp.e);

                            pos_t lce = lce_l + lce_r - 1;
                            qc_right = qc_right_tmp;

                            if (lce > f.len) {
                                f.len = lce;
                                f.src = C[idx_C.ssa(sample_point.y)] - lce_l + 1;

                                #ifndef NDEBUG
                                assert(f.src < i);

                                for (pos_t x = 0; x < f.len; x++) {
                                    assert(T[i + x] == T[f.src + x]);
                                }
                                #endif
                            }

                            return true;
                        }

                        if (log) {
                            num_range_queries++;
                            spa_range_sum += qc_left.e - qc_left.b + 1;
                            ssa_range_sum += qc_right_tmp.e - qc_right_tmp.b + 1;
                            time_range_queries += time_diff_ns(t3);
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