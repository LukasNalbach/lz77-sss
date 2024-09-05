#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
lz77_sss<pos_t>::factor lz77_sss<pos_t>::factorizer<tau>::longest_prev_occ(pos_t pos) {
    assert(gap_idx.pos() == pos);
    factor f{.src = char_to_uchar(T[pos]), .len = 0};

    for_constexpr<num_patt_lens - 1, - 1, - 1>([&](auto x) {
        if (f.len == 0) {
            pos_t src = gap_idx.template advance_and_get_occ<x>();

            if (src < pos && T[src] == T[pos]) {
                pos_t len = LCE_R(src, pos);

                #ifndef NDEBUG
                assert(src < pos);

                for (pos_t j = 0; j < len; j++) {
                    assert(T[pos + j] == T[src + j]);
                }
                #endif

                if (len > f.len) {
                    f.len = len;
                    f.src = src;
                }
            }
        } else {
            gap_idx.template advance<x>();
        }
    });

    gap_idx.inc_pos();
    return f;
}