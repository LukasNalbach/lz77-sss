#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
template <typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::factorize_greedy(out_it_t& out_it) {
    if (log) {
        std::cout << "factorizing" << std::flush;
    }
    
    pos_t p = 0;
    pos_t roll_threshold = 0;

    for_constexpr<0, num_patt_lens, 1>([&](auto j){
        roll_threshold += patt_lens[j];
    });

    roll_threshold /= num_patt_lens;
    
    for (pos_t i = 0; true;) {
        pos_t gap_end = LPF[p].beg;

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
                    if (i <= LPF[p].end) {
                        f.len -= i - gap_end;
                        i = gap_end;
                    } else {
                        do {
                            p++;
                        } while (LPF[p].end <= i);

                        while (gap_idx.pos() < gap_end) {
                            gap_idx.advance();
                        }

                        gap_end = LPF[p].beg;
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

                *out_it++ = f;
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
            .src = LPF[p].src + exc,
            .len = (LPF[p].end - LPF[p].beg) - exc
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

        *out_it++ = lpf;
        num_phr++;
        i += lpf.len;

        do {
            p++;
        } while (LPF[p].end <= i);
    }

    if (log) {
        time = log_runtime(time);
    }
}