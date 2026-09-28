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
void lz77_sss::factorizer<text_t>::exact_transformer::build_C()
{
    if (log) {
        std::cout << "setting delta = " << delta << std::endl;
        std::cout << "building C" << std::flush;
    }

    direct_ifstream aprx_ifile(aprx_file_name);
    std::istream_iterator<factor> aprx_it(aprx_ifile);
    aprx_it++;

    C.reset(num_aprx_fact + n / delta, n);
    C.push_back(0);

    num_par_sect = p == 1 ? 1 : (uint64_t{p} * par_sects_per_thr);
    par_sect.resize(num_par_sect + 1);
    par_sect[0] = sect_info_t {.beg = 0, .first_aprx = 0};
    par_sect[num_par_sect] = sect_info_t {.beg = n, .first_aprx = num_aprx_fact};

    uint64_t nxt_sect_aprx = num_aprx_fact / num_par_sect;
    uint64_t sect = 1;
    uint64_t cur_end = 0;
    uint64_t prev_end = 0;
    uint64_t prev_smpl;

    for (uint64_t k = 1; k < num_aprx_fact; k++) {
        prev_smpl = cur_end;
        factor f = *aprx_it++;
        cur_end += f.text_len();

        while (cur_end - prev_smpl > delta) {
            prev_smpl += delta;
            C.push_back(prev_smpl);
        }

        C.push_back(cur_end);

        if (k == nxt_sect_aprx) {
            par_sect[sect++] = sect_info_t {.beg = prev_end + 1, .first_aprx = k};
            nxt_sect_aprx = sect == num_par_sect ? num_aprx_fact : (sect * (num_aprx_fact / num_par_sect));
        }

        prev_end = cur_end;
    }

    C.finish();
    c = C.size();

    if (log) {
        record_phase_time("sample_set", time_diff_ns(time, now()));
        std::cout << " (" << format_size(C.size_in_bytes()) << ")";
        time = log_runtime(time);
        std::cout << "num. of samples / num. of aprx. factors = " << c / (double) num_aprx_fact << std::endl;
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::exact_transformer::build_idx_C()
{
    if (log) {
        std::cout << "building sample-index for C:" << std::endl;
    }

    uint64_t max_patt_len_left = delta;
    uint64_t max_smpl_len_right = get_max_smpl_len_right(n / (double) num_aprx_fact);
    const interval_samples mode = transf_mode == with_interval_samples ? interval_samples::use : interval_samples::skip;

    idx_C.build(T, n, C, LCE, mode, rks_sample_rate, p, log, max_patt_len_left, max_smpl_len_right);

    if (log) {
        std::cout << "size = " << format_size(idx_C.size_in_bytes());
        time = log_runtime(time);
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::exact_transformer::build_P()
{
    if (log) {
        std::cout << "building P" << std::flush;
    }

    P = bit_aligned_interleaved_vectors<3>(c, { c, c, kind.is_static() ? c : 0 });

    #pragma omp parallel for num_threads(p) schedule(dynamic, 65536)
    for (uint64_t i = 0; i < c; i++) {
        P.set_parallel<2>(i, i);
        P.set_parallel<0>(idx_C.pa(i), i);
        P.set_parallel<1>(idx_C.sa(i), i);
    }

    if (write_queries) {
        uint64_t num_points = c;
        result_log::queries_out.write((char*) &num_points, 8);

        for (uint64_t i = 0; i < c; i++) {
            uint64_t s_out = C[i];
            result_log::queries_out.write((char*) &s_out, 8);
        }

        for (uint64_t i = 0; i < c; i++) {
            struct pout_t { uint64_t x, y, weight; };
            pout_t p_out { .x = P.get<0>(i), .y = P.get<1>(i), .weight = P.get<2>(i) };
            result_log::queries_out.write((char*) &p_out, sizeof(pout_t));
        }
    }

    if (log) {
        record_phase_time("points", time_diff_ns(time, now()));
        std::cout << " (" << format_size(P.size_in_bytes()) << ")";
        time = log_runtime(time);
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::exact_transformer::build_Pi_Psi()
{
    if (log) {
        std::cout << "building Pi and Psi" << std::flush;
    }

    Pi = bit_aligned_vector(c, c);
    Psi = bit_aligned_vector(c, c);

    parallel_chunks(c, p, [&](uint64_t beg, uint64_t end) {
        Pi.fill(beg, end, [&](uint64_t i) { return P.get<1>(idx_C.pa(i)); });
        Psi.fill(beg, end, [&](uint64_t i) { return P.get<0>(idx_C.sa(i)); });
    });

    if (log) {
        record_phase_time("pi_psi", time_diff_ns(time, now()));
        std::cout << " (" << format_size(Pi.size_in_bytes() + Psi.size_in_bytes()) << ")";
        time = log_runtime(time);
    }
}

template <typename text_t>
inline void lz77_sss::factorizer<text_t>::exact_transformer::seek_x_c(uint64_t& x_c, uint64_t pos)
{
    while (x_c < c && C[x_c] < pos) {
        x_c++;
    }

    while (x_c > 0 && C[x_c - 1] >= pos) {
        x_c--;
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::exact_transformer::insert_points_before(uint64_t& x_r, uint64_t i)
{
    while (x_r < c && C[x_r] < i) {
        if (write_queries) {
            bool is_insert = true;
            char chr = T.char_at(C[x_r]);
            uint64_t weight = 0;
            uint64_t x_1 = P.get<0>(x_r);
            uint64_t x_2 = 0;
            uint64_t y_1 = P.get<1>(x_r);
            uint64_t y_2 = 0;

            result_log::queries_out.write((char*) &is_insert, 1);
            result_log::queries_out.write((char*) &chr, 1);
            result_log::queries_out.write((char*) &weight, 8);
            result_log::queries_out.write((char*) &x_1, 8);
            result_log::queries_out.write((char*) &x_2, 8);
            result_log::queries_out.write((char*) &y_1, 8);
            result_log::queries_out.write((char*) &y_2, 8);
        } else {
            R->insert(T[C[x_r]], point_t { .x = P.get<0>(x_r), .y = P.get<1>(x_r),
                .weight = P.get<2>(x_r) });
        }

        x_r++;
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::exact_transformer::find_close_sources(factor& f, uint64_t i, uint64_t e)
{
    uint64_t min_j = i <= delta ? 0 : (i - delta);

    for (uint64_t j = min_j; j < i; j++) {
        uint64_t lce = LCE_R(j, i);

        if (lce > f.len) [[unlikely]] {
            f.src = j;
            f.len = lce;
        }
    }

    if (f.len > e - i) [[unlikely]] {
        f.len = e - i;
    }
}

template <typename text_t>
bool lz77_sss::factorizer<text_t>::exact_transformer::intersect(
    const interval_t& pa_c_iv, const interval_t& sa_c_iv,
    uint64_t i, uint64_t j, uint64_t lce_l, uint64_t lce_r, uint64_t& x_c, factor& f)
{
    point_t point;
    bool found = false;

    uint64_t pa_c_width = pa_c_iv.e - pa_c_iv.b + 1;
    uint64_t sa_c_width = sa_c_iv.e - sa_c_iv.b + 1;

    if (std::min<uint64_t>(pa_c_width, sa_c_width) <= range_scan_threshold) {
        seek_x_c(x_c, j);

        if (pa_c_width <= sa_c_width) {
            for (uint64_t x = pa_c_iv.b; x <= pa_c_iv.e; x++) {
                if (idx_C.pa(x) < x_c && sa_c_iv.b <= Pi[x] && Pi[x] <= sa_c_iv.e) {
                    point.x = x;
                    point.y = Pi[x];
                    found = true;
                    break;
                }
            }
        } else {
            for (uint64_t y = sa_c_iv.b; y <= sa_c_iv.e; y++) {
                if (idx_C.sa(y) < x_c && pa_c_iv.b <= Psi[y] && Psi[y] <= pa_c_iv.e) {
                    point.x = Psi[y];
                    point.y = y;
                    found = true;
                    break;
                }
            }
        }
    } else if (kind.is_static()) {
        seek_x_c(x_c, j);

        if (write_queries) {
            bool is_insert = false;
            char chr = T.char_at(j);
            uint64_t weight = x_c;
            uint64_t x_1 = pa_c_iv.b;
            uint64_t x_2 = pa_c_iv.e;
            uint64_t y_1 = sa_c_iv.b;
            uint64_t y_2 = sa_c_iv.e;

            result_log::queries_out.write((char*) &is_insert, 1);
            result_log::queries_out.write((char*) &chr, 1);
            result_log::queries_out.write((char*) &weight, 8);
            result_log::queries_out.write((char*) &x_1, 8);
            result_log::queries_out.write((char*) &x_2, 8);
            result_log::queries_out.write((char*) &y_1, 8);
            result_log::queries_out.write((char*) &y_2, 8);
        }

        std::tie(point, found) = R->lighter_point_in_range(
            T[j], x_c,
            pa_c_iv.b, pa_c_iv.e,
            sa_c_iv.b, sa_c_iv.e);
    } else {
        std::tie(point, found) = R->point_in_range(
            T[j],
            pa_c_iv.b, pa_c_iv.e,
            sa_c_iv.b, sa_c_iv.e);
    }

    if (found) {
        #ifndef NDEBUG
        assert(pa_c_iv.b <= point.x && point.x <= pa_c_iv.e);
        assert(sa_c_iv.b <= point.y && point.y <= sa_c_iv.e);
        #endif

        uint64_t lce = lce_l + lce_r - 1;

        if (lce > f.len) {
            f.len = lce;
            f.src = C[idx_C.sa(point.y)] - lce_l + 1;

            #ifndef NDEBUG
            assert(f.src < i);

            for (uint64_t x = 0; x < f.len; x++) {
                assert(T[i + x] == T[f.src + x]);
            }
            #endif
        }
    }

    return found;
}

template <typename text_t>
lz77_sss::factor lz77_sss::factorizer<text_t>::exact_transformer::aprx_factor_at(uint64_t i)
{
    uint64_t sect = bin_search_max_leq<uint64_t, uint64_t>(
        i, 0, num_par_sect - 1, [&](uint64_t x) { return par_sect[x].beg; });
    uint64_t pos = par_sect[sect].beg;
    direct_ifstream aprx_ifile(aprx_file_name, 64 * 1024);
    aprx_ifile.seekg(par_sect[sect].first_aprx * factor::size_of(), std::ios::beg);
    std::istream_iterator<factor> aprx_it(aprx_ifile);
    factor f = *aprx_it++;

    while (pos + f.text_len() <= i) {
        pos += f.text_len();
        f = *aprx_it++;
    }

    if (!f.is_literal()) {
        uint64_t cut_left = i - pos;
        f.len = f.len - cut_left;
        f.src = f.src + cut_left;
    }

    return f;
}

template <typename text_t>
template <typename factor_at_t>
void lz77_sss::factorizer<text_t>::exact_transformer::
    combine_factorizations(factor_sink& output, factor_at_t factor_at)
{
    uint64_t pos = 0;
    num_fact = 0;

    auto emit = [&](factor f) {
        output(f);
        pos += f.text_len();
        num_fact++;
    };

    for (uint64_t sect = 0; sect < num_par_sect; sect++) {
        std::string fact_file_name_sect = fact_file_name + "_" + std::to_string(sect);
        const uint64_t e = par_sect[sect + 1].beg;
        const bool cut = sect + 1 < num_par_sect;
        uint64_t q = par_sect[sect].beg;

        {
            direct_ifstream fact_ifile(fact_file_name_sect);
            std::istream_iterator<factor> fact_it(fact_ifile);

            while (fact_it != std::istream_iterator<factor>()) {
                factor f = *fact_it++;
                const uint64_t beg = q;
                q += f.text_len();

                while (pos < beg) emit(factor_at(pos));
                if (beg < pos) continue;
                if (cut && q == e && !f.is_literal()) f = factor_at(pos);
                emit(f);
            }
        }

        std::filesystem::remove(fact_file_name_sect);
    }

    while (pos < n) emit(factor_at(pos));
}