#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::factorize_greedy_naive(std::function<void(factor)> out_it, std::function<lpf()> lpf_it) {
    if (log) {
        std::cout << "factorizing" << std::flush;
    }
    
    lpf p = lpf_it();
    lpf p_nxt;
    
    for (pos_t i = 0; true;) {
        pos_t gap_end = p.beg;

        if (i < gap_end) {
            if (gap_idx.pos() < i) {
                while (gap_idx.pos() < i) {
                    gap_idx.advance();
                }
            }

            do {
                factor f = longest_prev_occ(i);
                f.len = std::min<pos_t>(f.len, gap_end - i);
                out_it(f);
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
            .src = p.src,
            .len = p.end - p.beg
        };

        p_nxt = lpf_it();
        
        if (gap_idx.pos() == gap_end) {
            factor f = longest_prev_occ(gap_end);
            f.len = std::min<pos_t>(f.len, p_nxt.beg - gap_end);

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

        out_it(lpf);
        num_phr++;
        i += lpf.len;
        p = p_nxt;
    }

    if (log) {
        time = log_runtime(time);
    }
}