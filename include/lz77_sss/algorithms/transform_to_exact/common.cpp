#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::build_c(std::istream_iterator<factor>& ifile_approx_it) {
    if (log) {
        std::cout << "setting delta = " << delta << std::endl;
        std::cout << "building C" << std::flush;
    }

    C.reserve(num_phr + n / delta);
    C.emplace_back(0);
    pos_t end_lst = 0;
    pos_t end_nxt;

    for (pos_t i = 1; i < num_phr; i++) {
        end_nxt = end_lst + std::max<pos_t>(1,
            ifile_approx_it++->len);

        while (end_nxt - end_lst > delta) {
            end_lst += delta;
            C.emplace_back(end_lst);
        }

        C.emplace_back(end_nxt);
        end_lst = end_nxt;
    }

    C.shrink_to_fit();
    c = C.size();

    if (log) {
        std::cout << " (" << format_size(c * sizeof(pos_t)) << ")";
        time = log_runtime(time);
        std::cout << "c / num. of phrases = " << c / (double) num_phr << std::endl;
    }
}

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::build_idx_C() {
    if (log) {
        std::cout << "building sample-index for C:" << std::endl;
    }

    idx_C.build(T, C, LCE, delta, n / (1.0 * num_phr), transf_mode == optimized_with_samples, log);

    if (log) {
        std::cout << "size: " << format_size(idx_C.size_in_bytes());
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::build_p() {
    if (log) {
        std::cout << "building P" << std::flush;
    }

    P.reserve(c);

    for (sidx_t i = 0; i < c; i++) {
        if constexpr (is_static<range_ds_t>()) {
            P.emplace_back(point_t{.weight = i});
        } else {
            P.emplace_back(point_t{});
        }
    }

    for (sidx_t i = 0; i < c; i++) {
        P[idx_C.spa(i)].x = i;
        P[idx_C.ssa(i)].y = i;
    }

    if (log) {
        std::cout << " (" << format_size(sizeof(point_t) * c) << ")";
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::build_ps_sp() {
    if (log) {
        std::cout << "building PS and SP" << std::flush;
    }
    
    PS.reserve(c);

    for (sidx_t i = 0; i < c; i++) {
        PS.push_back(P[idx_C.spa(i)].y);
    }

    SP.reserve(c);

    for (sidx_t i = 0; i < c; i++) {
        SP.push_back(P[idx_C.ssa(i)].x);
    }

    if (log) {
        std::cout << " (" << format_size(2 * c * sizeof(sidx_t)) << ")";
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::adjust_xc(sidx_t& idx, pos_t pos) {
    while (idx < c && C[idx] < pos) {
        idx++;
    }

    while (idx > 0 && C[idx - 1] >= pos) {
        idx--;
    }
}

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::insert_points(sidx_t& x_c, pos_t i) {
    while (x_c < c && C[x_c] < i) {
        time_point_t t0;
        if (log) t0 = now();
        R.insert(P[x_c]);
        if (log) time_insert_points += time_diff_ns(t0);
        x_c++;
    }
}

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::handle_close_sources(factor& f, pos_t i) {
    pos_t min_j = i <= delta ? 0 : (i - delta);
    time_point_t t5;
    if (log) t5 = now();

    for (pos_t j = min_j; j < i; j++) {
        pos_t lce = LCE_R(j, i);
        
        if (lce > f.len) [[unlikely]] {
            f.src = j;
            f.len = lce;
        }
    }
    
    if (log) time_close_sources += time_diff_ns(t5);
}

template <typename pos_t>
template <uint64_t tau>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t, typename out_it_t>
bool lz77_sss<pos_t>::factorizer<tau>::exact_factorizer<sidx_t, transf_mode, range_ds_t, out_it_t>::intersect(
    const sxa_interval_t& spa_iv, const sxa_interval_t& ssa_iv,
    pos_t i, pos_t j, pos_t lce_l, pos_t lce_r, sidx_t& x_c, factor& f
) {
    time_point_t t3;
    point_t p;
    bool result = false;
    if (log) t3 = now();

    pos_t spa_rng = spa_iv.e - spa_iv.b + 1;
    pos_t ssa_rng = ssa_iv.e - ssa_iv.b + 1;

    if (transf_mode != naive && std::min<pos_t>(spa_rng, ssa_rng) <= range_scan_threshold) {
        adjust_xc(x_c, j);

        if (spa_rng <= ssa_rng) {
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
        } else {
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
        }
    } else if constexpr (is_static<range_ds_t>()) {
        adjust_xc(x_c, j);

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
        assert(spa_iv.b <= p.x && p.x <= spa_iv.e);
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