#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::greedy_phrase_selection(std::vector<lpf>& P) {
    if (P.empty()) {
        return;
    }

    ips4o::sort(P.begin(), P.end(), [](auto p1, auto p2) {
        return p1.beg < p2.beg ||
              (p1.beg == p2.beg && p1.end > p2.end);
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
void lz77_sss<pos_t>::factorizer<tau>::greedy_phrase_selection_external() {
    sel_lpf_file_name = "lpf_" + random_alphanumeric_string(10);
    std::ofstream lpf_ofile(sel_lpf_file_name);

    if (std::filesystem::file_size(lpf_file_name) == 0) {
        std::filesystem::remove(lpf_file_name);
        return;
    }

    std::ifstream lpf_ifile(lpf_file_name);
    std::istream_iterator<lpf> lpf_it_i(lpf_ifile);
    std::istream_iterator<lpf> lpf_it_x;
    std::ostream_iterator<lpf> lpf_out_it(lpf_ofile);

    pos_t k = 0;
    pos_t i = 1;
    pos_t p = std::filesystem::file_size(lpf_file_name) / sizeof(lpf);

    lpf lpf_k = *lpf_it_i++;
    lpf lpf_i = *lpf_it_i++;

    while (i < p && lpf_i.end < lpf_k.end) {
        i++;
        lpf_i = *lpf_it_i++;
    }

    while (i < p) {
        pos_t x = p;
        lpf_it_x = lpf_it_i;
        lpf lpf_x = *lpf_it_x++;

        if (i + 1 < p) {
            x = i + 1;
            
            while (x < p && lpf_x.beg <= lpf_k.end) {
                if (lpf_x.end > lpf_i.end) {
                    i = x;
                    lpf_i = lpf_x;
                    lpf_it_i = lpf_it_x;
                }

                x++;
                lpf_x = *lpf_it_x++;
            }

            if (lpf_i.end <= lpf_k.end) {
                i = x;
                lpf_i = lpf_x;
                lpf_it_i = lpf_it_x;
            }
        }

        if (i == p) {
            break;
        }

        if (lpf_i.beg < lpf_k.end) {
            lpf_k.end = lpf_i.beg;
        }

        if (lpf_k.end > lpf_k.beg) {
            k++;
            *lpf_out_it++ = lpf_k;
        }
        
        lpf_k = lpf_i;
        i = x;
        lpf_i = lpf_x;
        lpf_it_i = lpf_it_x;
    }

    std::filesystem::remove(lpf_file_name);
}

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::get_phrase_info() {
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

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::get_phrase_info_external() {
    num_lpf = std::filesystem::file_size(sel_lpf_file_name) / sizeof(lpf);

    if (num_lpf > 0) {
        std::ifstream lpf_ifile(sel_lpf_file_name);
        std::istream_iterator<lpf> lpf_it(lpf_ifile);
        lpf lpf_im1 = *lpf_it++;
        lpf lpf_i;
        
        if (num_lpf > 0) lpf_i = *lpf_it++;

        len_lpf_phr = lpf_im1.end - lpf_im1.beg;
        if (lpf_im1.beg > 0) num_gaps = 1;

        for (uint32_t i = 1; i < num_lpf; i++) {
            len_lpf_phr += lpf_i.end - lpf_i.beg;
            if (lpf_i.beg > lpf_im1.end) num_gaps++;
            lpf_im1 = lpf_i;
            lpf_i = *lpf_it++;
        }

        if (lpf_im1.end < n) num_gaps++;
    } else {
        num_gaps = 1;
    }

    len_gaps = n - len_lpf_phr;
}