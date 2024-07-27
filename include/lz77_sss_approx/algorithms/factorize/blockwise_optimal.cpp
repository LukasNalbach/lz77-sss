#pragma once

#include <lz77_sss_approx/lz77_sss_approx.hpp>

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
template <pos_t... patt_lens>
void lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::factorize_blockwise_optimal() {
    if (log) {
        std::cout << "initializing compact index" << std::flush;
    }

    rolling_hash_index_107<pos_t, patt_lens...> idx(T, target_index_size);

    if (log) {
        std::cout << " (size: " << format_size(idx.size_in_bytes()) << ")";
        time = log_runtime(time);
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

        assert(idx.pos() == block_start);

        for (pos_t i = block_start; i < block_end; i++) {
            for_constexpr<0, sizeof...(patt_lens), 1>([&](auto patt_len_idx) {
                constexpr pos_t len = ith_val<patt_len_idx>(patt_lens...);

                pos_t src = idx.template advance_and_get_occ<patt_len_idx>();

                if (src < i) {
                    pos_t lce_l = src == 0 ? 0 : LCE_L(src - 1, i - 1, i - block_start);
                    pos_t end = std::min<pos_t>(i + LCE_R(src, i), block_end);
                    pos_t beg = i - lce_l;

                    if (end - beg > 1) {
                        src -= lce_l;

                        #ifndef NDEBUG
                        assert(src < beg);
                        
                        for (pos_t j=0; j<end-beg; j++) {
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

            idx.inc_pos();
        }

        assert(idx.pos() == block_end);
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
            
            for (pos_t j=0; j<phr.end-phr.beg; j++) {
                assert(T[phr.src+j] == T[phr.beg+j]);
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

    assert(idx.pos() == n);

    if (log) {
        time = log_runtime(time);
        #ifndef NDEBUG
        std::cout << "rate of initialized values in the compact index: " << idx.rate_init() << std::endl;
        #endif
    }
}