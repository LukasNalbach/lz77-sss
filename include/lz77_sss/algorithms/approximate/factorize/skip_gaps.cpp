#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
template <typename lpf_it_t>
void lz77_sss<pos_t>::factorizer<tau>::factorize_skip_gaps(
    output_it_t& output,
    std::function<lpf_it_t()>& lpf_beg,
    std::function<lpf(lpf_it_t&)>& next_lpf)
{
    if (log) {
        std::cout << "outputting gapped factorization" << std::flush;
    }

    lpf_it_t lpf_it = lpf_beg();
    lpf phr;
    lpf phr_nxt = next_lpf(lpf_it);
    output({ .src = phr_nxt.beg, .len = 0 });

    while (phr_nxt.beg < n) {
        phr = phr_nxt;
        phr_nxt = next_lpf(lpf_it);
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