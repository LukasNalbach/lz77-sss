#pragma once

#include <lz77_sss/lz77_sss.hpp>
#include <gtl/btree.hpp>

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::build_LPF_all_external() {
    struct lpf_comp {
        bool operator()(const lpf &p1, const lpf &p2) const {
            return p1.beg < p2.beg ||
                (p1.beg == p2.beg && p1.end > p2.end);
        }
    };

    std::vector<gtl::btree_set<lpf, lpf_comp>> T_LPF(p);
    std::vector<std::vector<lpf>> LPF_buf(p);
    std::vector<std::ofstream> lpf_ofiles;
    lpf_ofiles.resize(p);
    pos_t tau_ = tau;

    for (uint16_t i_p = 0; i_p < p; i_p++) {
        open_lpf_ofile(lpf_ofiles[i_p], i_p);
    }

    build_LPF_all([&](uint16_t i_p, pos_t Si, lpf&& p){
        T_LPF[i_p].emplace_hint(T_LPF[i_p].end(), p);

        if (T_LPF[i_p].size() > max_buffered_lpf_phrases) {
            auto it = T_LPF[i_p].upper_bound(lpf {
                .beg = Si - tau_,
                .end = n,
                .src = 0
            });

            if (it != T_LPF[i_p].end()) {
                LPF_buf[i_p].insert(LPF_buf[i_p].end(),
                    T_LPF[i_p].begin(), it);

                T_LPF[i_p].erase(T_LPF[i_p].begin(), it);

                lpf_ofiles[i_p].write((char*) &LPF_buf[i_p][0],
                    LPF_buf[i_p].size() * sizeof(lpf));

                LPF_buf[i_p].clear();
            }
        }
    });

    #pragma omp parallel for num_threads(p)
    for (uint16_t i_p = 0; i_p < p; i_p++) {
        LPF_buf[i_p].insert(LPF_buf[i_p].end(),
        T_LPF[i_p].begin(), T_LPF[i_p].end());

        lpf_ofiles[i_p].write((char*) &LPF_buf[i_p][0],
            LPF_buf[i_p].size() * sizeof(lpf));

        LPF_buf[i_p].clear();
        lpf_ofiles[i_p].close();
    }
}

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::greedy_phrase_selection_external() {
    uint16_t i_p = omp_get_thread_num();
    std::ofstream lpf_ofile;
    open_sel_lpf_ofile(lpf_ofile, i_p);

    if (lpf_file_size(i_p) == 0) {
        remove_lpf_file(i_p);
        return;
    }

    std::ifstream lpf_ifile;
    std::istream_iterator<lpf> lpf_it_i, lpf_it_x;
    std::ostream_iterator<lpf> lpf_out_it;
    open_lpf_ifile(lpf_ifile, i_p);
    set_lpf_iterator(lpf_ifile, lpf_it_i);
    set_lpf_iterator(lpf_ofile, lpf_out_it);

    pos_t k = 0;
    pos_t i = 1;
    pos_t p = lpf_file_size(i_p);

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

    *lpf_out_it++ = lpf_k;

    lpf_ifile.close();
    lpf_ofile.close();
    remove_lpf_file(i_p);
}

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::get_phrase_info_external() {
    uint16_t i_p = omp_get_thread_num();

    pos_t num_gaps_thr = 0;
    pos_t len_lpf_phr_thr = 0;
    pos_t num_lpf_thr = 0;
    pos_t b = i_p * (n / p);
    pos_t e = i_p == p - 1 ? n : ((i_p + 1) * (n / p));
    pos_t i = 0;

    if (sel_lpf_file_size(i_p) != 0) {
        e = std::max<pos_t>(e, sel_lpf_file_back(i_p).end);
    }

    for (uint16_t j = 0; j < i_p; j++) {
        if (sel_lpf_file_size(j) != 0) {
            b = std::max<pos_t>(b, sel_lpf_file_back(j).end);
        }
    }

    pos_t lpf_ip_size = sel_lpf_file_size(i_p);
    std::ifstream lpf_ifile;
    std::istream_iterator<lpf> lpf_it;
    open_sel_lpf_ifile(lpf_ifile, i_p);
    set_lpf_iterator(lpf_ifile, lpf_it);
    lpf phr = *lpf_it++;
    lpf phr_lst;

    if (lpf_ip_size != 0) {
        while (i < lpf_ip_size && phr.end <= b) {
            phr = *lpf_it++;
            i++;
        }
    }

    num_lpf_thr = lpf_ip_size - i;

    if (num_lpf_thr > 0) {
        len_lpf_phr_thr += phr.end - std::max<pos_t>(phr.beg, b);
        if (phr.beg > b) num_gaps_thr = 1;
        phr_lst = phr;
        phr = *lpf_it++;
        i++;

        while (i < lpf_ip_size) {
            len_lpf_phr_thr += phr.end - phr.beg;
            if (phr.beg > phr_lst.end) num_gaps_thr++;
            phr_lst = phr;
            phr = *lpf_it++;
            i++;
        }

        if (phr.end < e) num_gaps_thr++;
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