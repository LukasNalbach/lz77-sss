#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
template <typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::factorize_blockwise_all(out_it_t& out_it, std::function<lpf()> lpf_it) {
    if (log) {
        std::cout << "factorizing" << std::flush;
    }

    pos_t blk_size = std::min<pos_t>(n, std::max<pos_t>(1000, n / 2000));
    pos_t blk_beg = 0;
    pos_t blk_end = blk_size;
    std::vector<lpf> P;
    lpf p = lpf_it();

    while (blk_beg < n) {
        while (p.beg < blk_end) {
            P.emplace_back(p);
            p = lpf_it();
        }

        if (!P.empty() && P.back().end > blk_end) {
            blk_end = P.back().end;
        }

        assert(gap_idx.pos() == blk_beg);

        for (pos_t i = blk_beg; i < blk_end; i++) {
            for_constexpr<0, num_patt_lens, 1>([&](auto x) {
                pos_t src = gap_idx.template advance_and_get_occ<x>();

                if (src < i) {
                    pos_t lce_l = src == 0 ? 0 : LCE_L(src - 1, i - 1, i - blk_beg);
                    pos_t end = std::min<pos_t>(i + LCE_R(src, i), blk_end);
                    pos_t beg = i - lce_l;

                    if (end - beg > 1) {
                        src -= lce_l;

                        #ifndef NDEBUG
                        assert(src < beg);
                        
                        for (pos_t j = 0; j < end - beg; j++) {
                            assert(T[src + j] == T[beg + j]);
                        }
                        #endif
                        
                        P.emplace_back(lpf {
                            .beg = beg,
                            .end = end,
                            .src = src
                        });
                    }
                }
            });

            gap_idx.inc_pos();
        }

        assert(gap_idx.pos() == blk_end);
        greedy_phrase_selection(P);
        pos_t pos = blk_beg;

        for (lpf phr : P) {
            while (pos < phr.beg) {
                *out_it++ = factor{char_to_uchar(T[pos]),0};
                num_phr++;
                pos++;
            }

            *out_it++ = factor{
                .src = phr.src,
                .len = phr.end - phr.beg
            };

            #ifndef NDEBUG
            assert(phr.src < pos);
            
            for (pos_t j = 0; j < phr.end - phr.beg; j++) {
                assert(T[phr.src + j] == T[phr.beg + j]);
            }
            #endif
            
            num_phr++;
            pos = phr.end;
        }

        while (pos < blk_end) {
            *out_it++ = factor{char_to_uchar(T[pos]),0};
            num_phr++;
            pos++;
        }

        P.clear();
        blk_beg = blk_end;
        blk_end = std::min<pos_t>(blk_beg + blk_size, n);
    }

    assert(gap_idx.pos() == n);

    if (log) {
        time = log_runtime(time);
    }
}