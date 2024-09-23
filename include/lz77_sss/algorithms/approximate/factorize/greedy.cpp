#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
template <typename lpf_it_t>
void lz77_sss<pos_t>::factorizer<tau>::factorize_greedy(
    output_it_t& output,
    std::function<lpf_it_t()>& lpf_beg,
    std::function<lpf(lpf_it_t&)>& next_lpf
) {
    if (log) {
        std::cout << "factorizing" << std::flush;
    }
    
    lpf_it_t lpf_it = lpf_beg();
    lpf p = next_lpf(lpf_it);
    
    for (pos_t i = 0; true;) {
        pos_t gap_end = p.beg;

        if (i < gap_end) {
            if (gap_idx.pos() < i) {
                if (i - gap_idx.pos() <= roll_threshold) {
                    do {
                        gap_idx.roll();
                    } while (gap_idx.pos() < i);
                } else {
                    gap_idx.reinit(i);
                }
            }

            do {
                factor f = longest_prev_occ(i);
                i += std::max<pos_t>(1, f.len);

                if (i > gap_end) {
                    if (i <= p.end) {
                        f.len -= i - gap_end;
                        i = gap_end;
                    } else {
                        do {
                            p = next_lpf(lpf_it);
                        } while (p.end <= i);

                        while (gap_idx.pos() < gap_end) {
                            gap_idx.advance();
                        }

                        gap_end = p.beg;
                    }
                }

                #ifndef NDEBUG
                pos_t i_ = i - std::max<pos_t>(1, f.len);

                assert(f.len == 0 ?
                    f.src == char_to_uchar(T[i_]) :
                    f.src < i_);

                for (pos_t j = 0; j < f.len; j++) {
                    assert(T[i_ + j] == T[f.src + j]);
                }
                #endif

                output(f);
                num_phr++;

                while (gap_idx.pos() < i) {
                    gap_idx.advance();
                }
            } while (i < gap_end);
        }

        if (i == n) {
            break;
        }

        pos_t exc = i - gap_end;

        factor lpf {
            .src = p.src + exc,
            .len = (p.end - p.beg) - exc
        };
        
        if (gap_idx.pos() == i) {
            factor f = longest_prev_occ(i);

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

        output(lpf);
        num_phr++;
        i += lpf.len;

        while (p.end <= i) {
            p = next_lpf(lpf_it);
        }
    }

    if (log) {
        time = log_runtime(time);
    }
}