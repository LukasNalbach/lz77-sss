#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <bool first_block>
inline lz77_sss<pos_t>::factor lz77_sss<pos_t>::factorizer<tau, char_t>::longest_prev_occ_par(fp_arr_t& fps, pos_t pos, pos_t end)
{
    factor f { .src = char_to_uchar(T[pos]), .len = 0 };

    for_constexpr<num_patt_lens - 1, -1, -1>([&](auto x) {
        if (f.len == 0) {
            pos_t src = par_gap_idx.template advance_and_get_occ<x, first_block>(fps, pos);

            if (src < pos && T[src] == T[pos]) {
                f.len = LCE_R(src, pos);
                f.src = src;

                #ifndef NDEBUG
                assert(f.src < pos);

                for (pos_t j = 0; j < f.len; j++) {
                    assert(T[pos + j] == T[f.src + j]);
                }
                #endif
            }
        } else {
            par_gap_idx.template advance<x>(fps, pos);
        }
    });

    pos_t len_max = end - pos;
    assert(len_max > 0);

    if (f.len > len_max) [[unlikely]] {
        f.len = len_max;
    }

    return f;
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <bool first_block, typename next_lpf_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::factorize_block(
    next_lpf_t next_lpf, pos_t blk_beg, pos_t blk_end)
{
    uint16_t i_p = omp_get_thread_num();

    block_info_t& info_beg = blk_info[blk_beg];
    block_info_t& info_end = blk_info[blk_end];
    pos_t beg = info_beg.beg;
    pos_t end = info_end.beg;
    fp_arr_t fps;
    lpf_pos_t lpf_pos;
    lpf_pos.i_p = info_beg.i_p;
    lpf_pos.i = info_beg.i;

    lpf phr = next_lpf(lpf_pos);
    pos_t pos_idx = beg;
    par_gap_idx.reinit(fps, beg);

    for (pos_t i = beg; true;) {
        pos_t gap_end = std::min<pos_t>(phr.beg, end);

        if (i < gap_end) {
            if (pos_idx < i) {
                if (i - pos_idx <= roll_threshold) {
                    do {
                        par_gap_idx.roll(fps, pos_idx);
                        pos_idx++;
                    } while (pos_idx < i);
                } else {
                    par_gap_idx.reinit(fps, i);
                    pos_idx = i;
                }
            }

            do {
                factor f = longest_prev_occ_par<first_block>(fps, i, end);
                pos_idx++;
                i += f.length();

                if (i > gap_end) {
                    if (i <= phr.end) {
                        f.len -= i - gap_end;
                        i = gap_end;
                    } else {
                        do {
                            phr = next_lpf(lpf_pos);
                        } while (phr.end <= i);

                        while (pos_idx < gap_end) {
                            par_gap_idx.advance(fps, pos_idx);
                            pos_idx++;
                        }

                        gap_end = std::min<pos_t>(phr.beg, end);
                    }
                }

                #ifndef NDEBUG
                pos_t i_ = i - f.length();

                assert(f.len == 0 ? f.src == char_to_uchar(T[i_]) : f.src < i_);

                for (pos_t j = 0; j < f.len; j++) {
                    assert(T[i_ + j] == T[f.src + j]);
                }
                #endif

                factors[i_p].emplace_back(f);

                while (pos_idx < i) {
                    par_gap_idx.advance(fps, pos_idx);
                    pos_idx++;
                }
            } while (i < gap_end);
        }

        if (i >= end) {
            break;
        }

        pos_t exc = i - gap_end;

        factor lpf {
            .src = phr.src + exc,
            .len = (phr.end - phr.beg) - exc
        };

        if (pos_idx == i) {
            factor f = longest_prev_occ_par<first_block>(fps, i, end);
            pos_idx++;

            if (f.len > lpf.len) {
                lpf = f;
            }
        }

        #ifndef NDEBUG
        assert(lpf.src < i);

        for (pos_t j = 0; j < lpf.len; j++) {
            assert(T[i + j] == T[lpf.src + j]);
        }
        #endif

        factors[i_p].emplace_back(lpf);
        i += lpf.len;
        assert(i <= end);

        while (phr.end <= i) {
            phr = next_lpf(lpf_pos);
            assert(phr.beg > 0);
        }
    }
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename output_fnc_t, typename lpf_beg_t, typename next_lpf_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::factorize_greedy_parallel(
    output_fnc_t output,
    lpf_beg_t lpf_beg,
    next_lpf_t next_lpf)
{
    if (log) {
        std::cout << "factorizing in parallel" << std::flush;
    }

    pos_t blk_size = std::max<pos_t>(min_par_gap_blk_size, (n / p) / max_par_gap_blks);
    pos_t num_blks = div_ceil(len_gaps, blk_size);
    blk_info.reserve(num_blks + 1);
    factors.resize(p);

    #ifndef NDEBUG
    pos_t cur_pos = 0;
    #endif

    {
        lpf_pos_t cur_lpf_pos = lpf_beg();
        lpf_pos_t lst_lpf_pos = cur_lpf_pos;
        lpf lpf_lst = next_lpf(cur_lpf_pos);
        lpf lpf_cur;
        pos_t cur_len_gaps = lpf_lst.beg;
        pos_t cur_blk_beg = 0;

        for (pos_t blk = 0; blk < num_blks; blk++) {
            while (cur_len_gaps <= cur_blk_beg) {
                lpf_cur = next_lpf(cur_lpf_pos);
                cur_len_gaps += lpf_cur.beg - lpf_lst.end;
                lst_lpf_pos = cur_lpf_pos;
                lpf_lst = lpf_cur;
            }

            blk_info.emplace_back(block_info_t {
                .i_p = lst_lpf_pos.i_p,
                .i = lst_lpf_pos.i,
                .beg = lpf_lst.beg - (cur_len_gaps - cur_blk_beg) });

            cur_blk_beg += blk_size;
        }

        blk_info.emplace_back(block_info_t { .beg = n });
    }

    for (pos_t cur_blk = 0; cur_blk < num_blks;) {
        pos_t blks = std::min<pos_t>(p, num_blks - cur_blk);
        par_gap_idx.overwrite(p);

        #pragma omp parallel num_threads(p)
        {
            uint16_t i_p = omp_get_thread_num();

            if (i_p < blks) {
                if (cur_blk == 0) {
                    if (i_p == 0) factorize_block<true, next_lpf_t>(next_lpf, 0, blks);
                } else {
                    pos_t blk_idx = cur_blk + i_p;
                    factorize_block<false, next_lpf_t>(next_lpf, blk_idx, blk_idx + 1);
                }
            }
        }

        for (auto& vec : factors) {
            for (factor f : vec) {
                output(f);

                #ifndef NDEBUG
                assert((f.len == 0 && (char) f.src == T[cur_pos]) || f.src < cur_pos);
                assert(f.len <= n - cur_pos);

                for (pos_t j = 0; j < f.len; j++) {
                    assert(T[f.src + j] == T[cur_pos + j]);
                }

                cur_pos += f.length();
                #endif
            }

            num_fact += vec.size();
            vec.clear();
        }

        cur_blk += blks;
    }

    blk_info.clear();
    blk_info.shrink_to_fit();
    factors.clear();
    factors.shrink_to_fit();

    if (log) {
        log_phase("factorize", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}