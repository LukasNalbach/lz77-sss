#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau, out_it_t>::greedy_phrase_selection(std::vector<lpf>& P) {
    ips4o::sort(P.begin(), P.end(), [](auto phr_i, auto phr_j) {
        return phr_i.beg < phr_j.beg ||
              (phr_i.beg == phr_j.beg && phr_i.end > phr_j.end);
    });

    pos_t k = 0;
    pos_t i = 1;
    pos_t p = P.size();

    while (i < p && P[i].end < P[k].end) i++;

    while (i < p) {
        pos_t x = p;

        if (i + 1 < p) {
            x = i + 1;
            
            while (x < p && P[x].beg <= P[k].end) {
                if (P[x].end > P[i].end) {
                    i = x;
                }

                x++;
            }

            if (P[i].end <= P[k].end) {
                i = x;
            }
        }

        if (i == p) {
            break;
        }

        if (P[i].beg < P[k].end) {
            P[k].end = P[i].beg;
        }

        if (P[k].end > P[k].beg) {
            k++;
        }
        
        P[k] = P[i];
        i = x;
    }
    
    P.resize(k + 1);

    #ifndef NDEBUG
    for (uint32_t i = 1; i < P.size(); i++) {
        assert(P[i - 1].end <= P[i].beg);
        assert(P[i - 1].beg < P[i].beg);
    }
    #endif
}

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau, out_it_t>::get_phrase_info() {
    num_lpf = LPF.size();

    if (num_lpf > 0) {
        len_lpf_phr = LPF[0].end - LPF[0].beg;
        if (LPF[0].beg > 0) num_gaps = 1;

        for (uint32_t i = 1; i < num_lpf; i++) {
            len_lpf_phr += LPF[i].end - LPF[i].beg;
            if (LPF[i].beg > LPF[i - 1].end) num_gaps++;
        }

        if (LPF[num_lpf - 1].end < n) num_gaps++;
    } else {
        num_gaps = 1;
    }

    len_gaps = n - len_lpf_phr;
}