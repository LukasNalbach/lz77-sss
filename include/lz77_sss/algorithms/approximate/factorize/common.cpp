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