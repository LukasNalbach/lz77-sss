#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::build_c()
{
    if (log) {
        std::cout << "setting delta = " << delta << std::endl;
        std::cout << "building C" << std::flush;
    }

    std::ifstream aprx_ifile(aprx_file_name);
    std::istream_iterator<factor> aprx_it(aprx_ifile);
    aprx_it++;

    C.reserve(num_fact + n / delta);
    C.emplace_back(0);

    par_sect.resize(p + 1);
    par_sect[0] = sect_info_t {.beg = 0, .phr_idx = 0};
    par_sect[p] = sect_info_t {.beg = n, .phr_idx = num_fact};

    sidx_t phr_nxt = num_fact / p;
    uint16_t i_p = 1;
    pos_t end_cur = 0;
    pos_t end_lst = 0;
    pos_t smpl_lst;

    for (sidx_t phr = 1; phr < num_fact; phr++) {
        smpl_lst = end_cur;
        factor fac = *aprx_it++;
        end_cur += fac.length();

        while (end_cur - smpl_lst > delta) {
            smpl_lst += delta;
            C.emplace_back(smpl_lst);
        }

        C.emplace_back(end_cur);

        if (phr == phr_nxt) {
            par_sect[i_p++] = sect_info_t {.beg = end_lst + 1, .phr_idx = phr};
            phr_nxt = i_p == p ? num_fact : (i_p * (num_fact / p));
        }

        end_lst = end_cur;
    }

    C.shrink_to_fit();
    c = C.size();

    if (log) {
        log_phase("sample_set", time_diff_ns(time, now()));
        std::cout << " (" << format_size(c * sizeof(pos_t)) << ")";
        time = log_runtime(time);
        std::cout << "num. of samples / num. of aprx. factors = " << c / (double) num_fact << std::endl;
    }
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::build_idx_C()
{
    if (log) {
        std::cout << "building sample-index for C:" << std::endl;
    }

    pos_t max_patt_len_left = delta;
    pos_t max_smpl_len_right = get_max_smpl_len_right(n / (double) num_fact);
    constexpr bool build_interval_samples = transf_mode == with_samples;

    idx_C.build(T, n, C, LCE, build_interval_samples, rks_sample_rate, p, log, max_patt_len_left, max_smpl_len_right);

    if (log) {
        std::cout << "size: " << format_size(idx_C.size_in_bytes());
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::build_p()
{
    if (log) {
        std::cout << "building P" << std::flush;
    }

    P.resize(c, point_t { });

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < c; i++) {
        if constexpr (range_ds_t<sidx_t>::is_static()) {
            P[i].weight = i;
        }

        P[idx_C.pa_s(i)].x = i;
        P[idx_C.sa_s(i)].y = i;
    }

    #if defined(GEN_RANGE_QUERIES)
    uint64_t num_points = c;
    queries_file.write((char*) &num_points, 8);

    for (pos_t s : C) {
        uint64_t s_out = s;
        queries_file.write((char*) &s_out, 8);
    }

    for (point_t p : P) {
        using pout_t = static_weighted_range<uint64_t>::point_t;
        pout_t p_out { .x = p.x, .y = p.y, .weight = p.weight };
        queries_file.write((char*) &p_out, sizeof(pout_t));
    }
    #endif

    if (log) {
        log_phase("points", time_diff_ns(time, now()));
        std::cout << " (" << format_size(sizeof(point_t) * c) << ")";
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::build_pi_psi()
{
    if (log) {
        std::cout << "building Pi and Psi" << std::flush;
    }

    no_init_resize(Pi, c);
    no_init_resize(Psi, c);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < c; i++) {
        Pi[i] = P[idx_C.pa_s(i)].y;
    }

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < c; i++) {
        Psi[i] = P[idx_C.sa_s(i)].x;
    }

    if (log) {
        log_phase("ypsilon", time_diff_ns(time, now()));
        std::cout << " (" << format_size(2 * c * sizeof(sidx_t)) << ")";
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
inline void lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::adjust_xc(sidx_t& idx, pos_t pos)
{
    while (idx < c && C[idx] < pos) {
        idx++;
    }

    while (idx > 0 && C[idx - 1] >= pos) {
        idx--;
    }
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::insert_points(sidx_t& x_r, pos_t i)
{
    while (x_r < c && C[x_r] < i) {
        #if not defined(GEN_RANGE_QUERIES)
        if constexpr (range_ds_t<sidx_t>::is_decomposed()) {
            R.insert(T[C[x_r]], P[x_r]);
        } else {
            R.insert(P[x_r]);
        }
        #endif

        #if defined(GEN_RANGE_QUERIES)
        bool is_insert = true;
        char c = T[C[x_r]];
        uint64_t weight = 0;
        uint64_t x_1 = P[x_r].x;
        uint64_t x_2 = 0;
        uint64_t y_1 = P[x_r].y;
        uint64_t y_2 = 0;

        queries_file.write((char*) &is_insert, 1);
        queries_file.write((char*) &c, 1);
        queries_file.write((char*) &weight, 8);
        queries_file.write((char*) &x_1, 8);
        queries_file.write((char*) &x_2, 8);
        queries_file.write((char*) &y_1, 8);
        queries_file.write((char*) &y_2, 8);
        #endif

        x_r++;
    }
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::find_close_sources(factor& f, pos_t i, pos_t e)
{
    pos_t min_j = i <= delta ? 0 : (i - delta);

    for (pos_t j = min_j; j < i; j++) {
        pos_t lce = LCE_R(j, i);

        if (lce > f.len) [[unlikely]] {
            f.src = j;
            f.len = lce;
        }
    }

    if (f.len > e - i) [[unlikely]] {
        f.len = e - i;
    }
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
bool lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::intersect(
    const interval_t& pa_c_iv, const interval_t& sa_c_iv,
    pos_t i, pos_t j, pos_t lce_l, pos_t lce_r, sidx_t& x_c, factor& f)
{
    point_t p;
    bool result = false;

    pos_t pa_s_rng = pa_c_iv.e - pa_c_iv.b + 1;
    pos_t sa_s_rng = sa_c_iv.e - sa_c_iv.b + 1;

    if (transf_mode != naive && std::min<pos_t>(pa_s_rng, sa_s_rng) <= range_scan_threshold) {
        adjust_xc(x_c, j);

        if (pa_s_rng <= sa_s_rng) {
            for (sidx_t x = pa_c_iv.b; x <= pa_c_iv.e; x++) {
                if (idx_C.pa_s(x) < x_c && sa_c_iv.b <= Pi[x] && Pi[x] <= sa_c_iv.e) {
                    p.x = x;
                    p.y = Pi[x];
                    result = true;
                    break;
                }
            }
        } else {
            for (sidx_t y = sa_c_iv.b; y <= sa_c_iv.e; y++) {
                if (idx_C.sa_s(y) < x_c && pa_c_iv.b <= Psi[y] && Psi[y] <= pa_c_iv.e) {
                    p.x = Psi[y];
                    p.y = y;
                    result = true;
                    break;
                }
            }
        }
    } else if constexpr (range_ds_t<sidx_t>::is_static()) {
        adjust_xc(x_c, j);

        #if defined(GEN_RANGE_QUERIES)
        bool is_insert = false;
        char c = T[j];
        uint64_t weight = x_c;
        uint64_t x_1 = pa_c_iv.b;
        uint64_t x_2 = pa_c_iv.e;
        uint64_t y_1 = sa_c_iv.b;
        uint64_t y_2 = sa_c_iv.e;

        queries_file.write((char*) &is_insert, 1);
        queries_file.write((char*) &c, 1);
        queries_file.write((char*) &weight, 8);
        queries_file.write((char*) &x_1, 8);
        queries_file.write((char*) &x_2, 8);
        queries_file.write((char*) &y_1, 8);
        queries_file.write((char*) &y_2, 8);
        #endif

        if constexpr (range_ds_t<sidx_t>::is_decomposed()) {
            std::tie(p, result) = R.lighter_point_in_range(
                T[j], x_c,
                pa_c_iv.b, pa_c_iv.e,
                sa_c_iv.b, sa_c_iv.e);
        } else {
            std::tie(p, result) = R.lighter_point_in_range(
                x_c,
                pa_c_iv.b, pa_c_iv.e,
                sa_c_iv.b, sa_c_iv.e);
        }
    } else {
        if constexpr (range_ds_t<sidx_t>::is_decomposed()) {
            std::tie(p, result) = R.point_in_range(
                T[j],
                pa_c_iv.b, pa_c_iv.e,
                sa_c_iv.b, sa_c_iv.e);
        } else {
            std::tie(p, result) = R.point_in_range(
                pa_c_iv.b, pa_c_iv.e,
                sa_c_iv.b, sa_c_iv.e);
        }
    }

    if (result) {
        #ifndef NDEBUG
        assert(pa_c_iv.b <= p.x && p.x <= pa_c_iv.e);
        assert(sa_c_iv.b <= p.y && p.y <= sa_c_iv.e);
        #endif

        pos_t lce = lce_l + lce_r - 1;

        if (lce > f.len) {
            f.len = lce;
            f.src = C[idx_C.sa_s(p.y)] - lce_l + 1;

            #ifndef NDEBUG
            assert(f.src < i);

            for (pos_t x = 0; x < f.len; x++) {
                assert(T[i + x] == T[f.src + x]);
            }
            #endif
        }
    }

    return result;
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
template <typename output_fnc_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::
    combine_factorizations(output_fnc_t output)
{
    for (uint16_t i_p = 0; i_p < p; i_p++) {
        std::string fact_file_name_thr = fact_file_name + "_" + std::to_string(i_p);
        std::ifstream fact_ifile(fact_file_name_thr);
        std::istream_iterator<factor> fact_it(fact_ifile);

        while (fact_it != std::istream_iterator<factor>()) {
            output(*fact_it++);
        }

        fact_ifile.close();
        std::filesystem::remove(fact_file_name_thr);
    }
}