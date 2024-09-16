#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::factorize_skip_gaps(std::function<void(factor)> out_it, std::function<lpf()> lpf_it) {
    if (log) {
        std::cout << "outputting gapped factorization" << std::flush;
    }

    lpf phr;
    lpf phr_nxt = lpf_it();
    out_it({.src = phr_nxt.beg, .len = 0});

    while (phr_nxt.beg < n) {
        phr = phr_nxt;
        phr_nxt = lpf_it();

        out_it({
            .src = phr.src,
            .len = phr.end - phr.beg
        });

        if (phr_nxt.beg > phr.end) {
            out_it({
                .src = phr_nxt.beg - phr.end,
                .len = 0
            });
        }
    }

    if (log) {
        time = log_runtime(time);
    }
}