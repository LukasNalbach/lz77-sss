#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename output_fnc_t, typename lpf_beg_t, typename next_lpf_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::factorize_greedy_naive(
    output_fnc_t output,
    lpf_beg_t lpf_beg,
    next_lpf_t next_lpf)
{
    if (log) {
        std::cout << "factorizing sequentially" << std::flush;
    }

    lpf_pos_t lpf_pos = lpf_beg();
    lpf p = next_lpf(lpf_pos);
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
                output(f);
                num_fact++;
                i += f.length();

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

        p_nxt = next_lpf(lpf_pos);

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

        output(lpf);
        num_fact++;
        i += lpf.len;
        p = p_nxt;
    }

    if (log) {
        log_phase("factorize", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}