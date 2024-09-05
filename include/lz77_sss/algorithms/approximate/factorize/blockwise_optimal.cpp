#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
template <typename out_it_t>
void lz77_sss<pos_t>::factorizer<tau>::factorize_blockwise_optimal(out_it_t& out_it) {
    if (log) {
        std::cout << "factorizing" << std::flush;
    }

    pos_t block_size = std::min<pos_t>(n, std::max<pos_t>(1000, n / 2000));
    pos_t block_start = 0;
    pos_t block_end = block_size;
    std::vector<lpf> phrases;
    pos_t cur_lpf = 0;

    while (block_start < n) {
        while (cur_lpf < LPF.size() && LPF[cur_lpf].beg < block_end) {
            phrases.emplace_back(LPF[cur_lpf]);
            cur_lpf++;
        }

        if (!phrases.empty() && phrases.back().end > block_end) {
            block_end = phrases.back().end;
        }

        assert(gap_idx.pos() == block_start);

        for (pos_t i = block_start; i < block_end; i++) {
            for_constexpr<0, num_patt_lens, 1>([&](auto x) {
                pos_t src = gap_idx.template advance_and_get_occ<x>();

                if (src < i) {
                    pos_t lce_l = src == 0 ? 0 : LCE_L(src - 1, i - 1, i - block_start);
                    pos_t end = std::min<pos_t>(i + LCE_R(src, i), block_end);
                    pos_t beg = i - lce_l;

                    if (end - beg > 1) {
                        src -= lce_l;

                        #ifndef NDEBUG
                        assert(src < beg);
                        
                        for (pos_t j = 0; j < end - beg; j++) {
                            assert(T[src + j] == T[beg + j]);
                        }
                        #endif
                        
                        phrases.emplace_back(lpf {
                            .beg = beg,
                            .end = end,
                            .src = src
                        });
                    }
                }
            });

            gap_idx.inc_pos();
        }

        assert(gap_idx.pos() == block_end);
        greedy_phrase_selection(phrases);
        pos_t pos = block_start;

        for (lpf phr : phrases) {
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

        while (pos < block_end) {
            *out_it++ = factor{char_to_uchar(T[pos]),0};
            num_phr++;
            pos++;
        }

        phrases.clear();
        block_start = block_end;
        block_end = std::min<pos_t>(block_start + block_size, n);
    }

    assert(gap_idx.pos() == n);

    if (log) {
        time = log_runtime(time);
    }
}