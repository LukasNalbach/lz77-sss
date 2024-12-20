#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
inline lz77_sss<pos_t>::factor lz77_sss<pos_t>::factorizer<tau>::longest_prev_occ(pos_t pos)
{
    assert(gap_idx.pos() == pos);
    factor f { .src = char_to_uchar(T[pos]), .len = 0 };

    for_constexpr<num_patt_lens - 1, -1, -1>([&](auto x) {
        if (f.len == 0) {
            pos_t src = gap_idx.template advance_and_get_occ<x>();

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
            gap_idx.template advance<x>();
        }
    });

    gap_idx.inc_pos();
    return f;
}

template <typename pos_t>
template <uint64_t tau>
template <factorize_mode fact_mode>
void lz77_sss<pos_t>::factorizer<tau>::factorize(output_it_t& output)
{
    std::function<lpf_arr_it_t()> lpf_beg = [&]() {
        lpf_arr_it_t lpf_it { .i_p = 0, .i = 0 };
        while (LPF[lpf_it.i_p].empty()) lpf_it.i_p++;
        return lpf_it;
    };

    std::function<lpf(lpf_arr_it_t&)> next_lpf = [&](lpf_arr_it_t& lpf_it) {
        uint16_t& i_p = lpf_it.i_p;
        uint32_t& i = lpf_it.i;

        lpf phr = LPF[i_p][i++];

        if (i == LPF[i_p].size()) [[unlikely]] {
            if (i_p < p - 1 && p > 1) [[likely]] {
                do {
                    i_p++;
                    i = 0;

                    while (i < LPF[i_p].size() && LPF[i_p][i].end <= phr.end) {
                        i++;
                    }

                    lpf& phr_i = LPF[i_p][i];

                    if (i < LPF[i_p].size() && phr_i.beg < phr.end) {
                        pos_t offs = phr.end - phr_i.beg;
                        phr_i.beg += offs;
                        phr_i.src += offs;
                    }
                } while (i == LPF[i_p].size() && i_p < p - 1);
            }
        }

        return phr;
    };

    if (fact_mode == greedy && par_gap_idx.size() != 0) {
        factorize_greedy_parallel(output, lpf_beg, next_lpf);
    } else {
        factorize_sequential<fact_mode>(output, lpf_beg, next_lpf);
    }
}