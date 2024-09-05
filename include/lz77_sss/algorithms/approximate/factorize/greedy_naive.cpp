#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
template <typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::factorize_greedy_naive(out_it_t& out_it) {
    if (log) {
        std::cout << "factorizing" << std::flush;
    }
    
    pos_t p = 0;
    pos_t gap_start = 0;
    pos_t patt_len_sum = 0;

    for_constexpr<0, num_patt_lens, 1>([&](auto i){
        patt_len_sum += patt_lens[i];
    });
    
    for (pos_t i = 0; true;) {
        pos_t gap_end = LPF[p].beg;

        if (gap_start < gap_end) {
            if (gap_idx.pos() < gap_start) {
                while (gap_idx.pos() < gap_start) {
                    gap_idx.advance();
                }
            }

            do {
                factor f = longest_prev_occ(gap_start);
                f.len = std::min<pos_t>(f.len, gap_end - gap_start);
                *out_it++ = f;
                num_phr++;
                gap_start += std::max<pos_t>(1, f.len);

                while (gap_idx.pos() < gap_start) {
                    gap_idx.advance();
                }
            } while (gap_start < gap_end);
        }

        if (gap_end == n) {
            break;
        }

        factor lpf {
            .src = LPF[p].src,
            .len = LPF[p].end - LPF[p].beg
        };
        
        if (gap_idx.pos() == gap_end) {
            pos_t beg_nxt_lpf = p + 1 == num_lpf ? n : LPF[p + 1].beg;
            factor f = longest_prev_occ(gap_end);
            f.len = std::min<pos_t>(f.len, beg_nxt_lpf - gap_end);

            if (f.len > lpf.len) {
                lpf = f;
            }
        }

        #ifndef NDEBUG
        assert(lpf.src < gap_end);

        for (pos_t j = 0; j < lpf.len; j++) {
            assert(T[gap_end + j] == T[lpf.src + j]);
        }
        #endif

        *out_it++ = lpf;
        num_phr++;
        gap_start += lpf.len;
        p++;
    }

    if (log) {
        time = log_runtime(time);
    }
}