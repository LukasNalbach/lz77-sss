#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::greedy_phrase_selection(std::vector<lpf>& P)
{
    if (P.empty()) {
        return;
    }

    ips4o::sort(P.begin(), P.end(), [](auto p1, auto p2) {
        return p1.beg < p2.beg || (p1.beg == p2.beg && p1.end > p2.end);
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

    for (uint32_t i = 0; i < P.size(); i++) {
        assert(P[i].beg < P[i].end);
    }
    #endif
}

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::get_phrase_info()
{
#pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t num_gaps_thr = 0;
        pos_t len_lpf_phr_thr = 0;
        pos_t num_lpf_thr = 0;
        pos_t b = i_p * (n / p);
        pos_t e = i_p == p - 1 ? n : ((i_p + 1) * (n / p));
        pos_t i = 0;

        if (!LPF[i_p].empty()) {
            e = std::max<pos_t>(e, LPF[i_p].back().end);
        }

        for (uint16_t j = 0; j < i_p; j++) {
            if (!LPF[j].empty()) {
                b = std::max<pos_t>(b, LPF[j].back().end);
            }
        }

        if (!LPF[i_p].empty()) {
            while (i < LPF[i_p].size() && LPF[i_p][i].end <= b) {
                i++;
            }
        }

        num_lpf_thr = LPF[i_p].size() - i;

        if (num_lpf_thr > 0) {
            len_lpf_phr_thr += LPF[i_p][i].end - std::max<pos_t>(LPF[i_p][i].beg, b);
            if (LPF[i_p][i].beg > b)
                num_gaps_thr = 1;
            i++;

            while (i < LPF[i_p].size()) {
                len_lpf_phr_thr += LPF[i_p][i].end - LPF[i_p][i].beg;
                if (LPF[i_p][i].beg > LPF[i_p][i - 1].end)
                    num_gaps_thr++;
                i++;
            }

            if (LPF[i_p].back().end < e)
                num_gaps_thr++;
        } else {
            num_gaps_thr = 1;
        }

        #pragma omp atomic
        num_gaps += num_gaps_thr;
        #pragma omp atomic
        num_lpf += num_lpf_thr;
        #pragma omp atomic
        len_lpf_phr += len_lpf_phr_thr;
    }
}