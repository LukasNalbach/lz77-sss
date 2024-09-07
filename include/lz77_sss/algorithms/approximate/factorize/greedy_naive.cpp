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
    
    for (pos_t i = 0; true;) {
        pos_t gap_end = LPF[p].beg;

        if (i < gap_end) {
            if (gap_idx.pos() < i) {
                while (gap_idx.pos() < i) {
                    gap_idx.advance();
                }
            }

            do {
                factor f = longest_prev_occ(i);
                f.len = std::min<pos_t>(f.len, gap_end - i);
                *out_it++ = f;
                num_phr++;
                i += std::max<pos_t>(1, f.len);

                while (gap_idx.pos() < i) {
                    gap_idx.advance();
                }
            } while (i < gap_end);
        }

        if (gap_end == n) {
            break;
        }

        factor lpf {
            .src = LPF[p].src,
            .len = LPF[p].end - LPF[p].beg
        };
        
        if (gap_idx.pos() == gap_end) {
            factor f = longest_prev_occ(gap_end);
            f.len = std::min<pos_t>(f.len, LPF[p + 1].beg - gap_end);

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
        i += lpf.len;
        p++;
    }

    if (log) {
        time = log_runtime(time);
    }
}