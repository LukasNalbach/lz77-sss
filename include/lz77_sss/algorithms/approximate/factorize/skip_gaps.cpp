#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename char_t>
template <typename output_fnc_t, typename lpf_beg_t, typename next_lpf_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::factorize_skip_gaps(
    output_fnc_t output,
    lpf_beg_t lpf_beg,
    next_lpf_t next_lpf)
{
    if (log) {
        std::cout << "outputting gapped factorization" << std::flush;
    }

    lpf_pos_t lpf_pos = lpf_beg();
    lpf phr;
    lpf phr_nxt = next_lpf(lpf_pos);
    output({ .src = phr_nxt.beg, .len = 0 });

    while (phr_nxt.beg < n) {
        phr = phr_nxt;
        phr_nxt = next_lpf(lpf_pos);
        output({ .src = phr.src, .len = phr.end - phr.beg });

        if (phr_nxt.beg > phr.end) {
            output({ .src = phr_nxt.beg - phr.end, .len = 0 });
        }
    }

    if (log) {
        log_phase("output_gapped", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}