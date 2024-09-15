#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
inline lz77_sss<pos_t>::factor lz77_sss<pos_t>::factorizer<tau>::longest_prev_occ(pos_t pos) {
    assert(gap_idx.pos() == pos);
    factor f{.src = char_to_uchar(T[pos]), .len = 0};

    for_constexpr<num_patt_lens - 1, - 1, - 1>([&](auto x) {
        if (f.len == 0) {
            pos_t src = gap_idx.template advance_and_get_occ<x>();

            if (src < pos && T[src] == T[pos]) {
                f.len = LCE_R(src, pos);
                f.src = src;

                #ifndef NDEBUG
                assert(f.src < pos);

                for (pos_t j = 0; j < f.len; j++) {
                    assert(T[pos + j] == T[f.src + j]);
                }
                #endif
            }
        } else {
            gap_idx.template advance<x>();
        }
    });

    gap_idx.inc_pos();
    return f;
}

template <typename pos_t>
template <uint64_t tau>
template <factorize_mode fact_mode, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::factorize_internal(out_it_t& out_it) {
    uint16_t i_p = 0;
    uint32_t i = 0;
    while(LPF[i_p].empty()) i_p++;

    factorize<fact_mode>(out_it, [&](){
        lpf phr = LPF[i_p][i++];
        
        if (i == LPF[i_p].size()) [[unlikely]] {
            if (i_p < p - 1 && p > 1) [[likely]] {
                do {
                    i_p++;
                    i = 0;

                    while (i < LPF[i_p].size() && LPF[i_p][i].end <= phr.end) {
                        i++;
                    }

                    if (i < LPF[i_p].size() && LPF[i_p][i].beg < phr.end) {
                        pos_t offs = phr.end - LPF[i_p][i].beg;
                        LPF[i_p][i].beg += offs;
                        LPF[i_p][i].src += offs;
                    }
                } while (i == LPF[i_p].size() && i_p < p - 1);
            }
        }

        return phr;
    });
}

template <typename pos_t>
template <uint64_t tau>
template <factorize_mode fact_mode, typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::factorize_external(out_it_t& out_it) {
    uint16_t i_p = 0;
    while (sel_lpf_file_size(i_p) == 0) i_p++;
    std::ifstream lpf_ifile;
    std::istream_iterator<lpf> lpf_it;
    open_sel_lpf_ifile(lpf_ifile, i_p);
    set_lpf_iterator(lpf_ifile, lpf_it);
    constexpr auto end_it = std::istream_iterator<lpf>();
    lpf phr_nxt = *lpf_it++;

    factorize<fact_mode>(out_it, [&](){
        lpf phr = phr_nxt;

        if (lpf_it == end_it) [[unlikely]] {
            if (i_p < p - 1 && p > 1) [[likely]] {
                do {
                    i_p++;
                    lpf_ifile.close();
                    open_sel_lpf_ifile(lpf_ifile, i_p);
                    set_lpf_iterator(lpf_ifile, lpf_it);

                    while (lpf_it != end_it) {
                        phr_nxt = *lpf_it++;

                        if (phr_nxt.end > phr.end) {
                            break;
                        }
                    }

                    if (lpf_it != end_it && phr_nxt.beg < phr.end) {
                        pos_t offs = phr.end - phr_nxt.beg;
                        phr_nxt.beg += offs;
                        phr_nxt.src += offs;
                    }
                } while (lpf_it == end_it && i_p < p - 1);
            }
        } else {
            phr_nxt = *lpf_it++;
        }

        return phr;
    });

    for (uint16_t i_p = 0; i_p < p; i_p++) {
        remove_sel_lpf_file(i_p);
    }
}