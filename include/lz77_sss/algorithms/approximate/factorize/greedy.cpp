#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
template <pos_t... patt_lens>
lz77_sss<pos_t>::factor lz77_sss<pos_t>::factorizer<tau, out_it_t>::longest_prev_occ(
    rolling_hash_index_107<pos_t, patt_lens...>& idx, pos_t pos, pos_t len_max
) {
    assert(idx.pos() == pos);
    factor f{0, 0};

    for_constexpr<sizeof...(patt_lens) - 1, - 1, - 1>([&](auto patt_len_idx) {
        if (f.len == 0) {
            pos_t src = idx.template advance_and_get_occ<patt_len_idx>();

            if (src < pos && T[src] == T[pos]) {
                pos_t len = std::min<pos_t>(LCE_R(src, pos), len_max);

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
            idx.template advance<patt_len_idx>();
        }
    });

    idx.inc_pos();

    if (f.len == 0) {
        f.src = char_to_uchar(T[pos]);
    }

    return f;
}

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
template <bool skip_phrases, pos_t... patt_lens>
void lz77_sss<pos_t>::factorizer<tau, out_it_t>::factorize_greedy(out_it_t& out_it) {
    if (log) {
        std::cout << "initializing rolling hash index" << std::flush;
    }

    rolling_hash_index_107<pos_t, patt_lens...> idx(T, target_index_size);

    if (log) {
        std::cout << " (size: " << format_size(idx.size_in_bytes()) << ")";
        time = log_runtime(time);
        std::cout << "factorizing" << std::flush;
    }
    
    pos_t cur_lpf = 0;
    pos_t gap_start = 0;
    
    while (true) {
        pos_t gap_end = cur_lpf == num_lpf ? n : LPF[cur_lpf].beg;

        if (gap_start < gap_end) {
            if (idx.pos() < gap_start) {
                if constexpr (skip_phrases) {
                    if ((gap_start - idx.pos()) * (2 * sizeof...(patt_lens)) <
                        constexpr_sum<pos_t>(patt_lens...)
                    ) {
                        while (idx.pos() < gap_start) {
                            idx.roll();
                        }
                    } else {
                        idx.reinit(gap_start);
                    }
                } else {
                    while (idx.pos() < gap_start) {
                        idx.advance();
                    }
                }
            }

            do {
                factor f = longest_prev_occ(idx, gap_start, gap_end - gap_start);
                *out_it++ = f;
                num_phr++;
                gap_start += std::max<pos_t>(1, f.len);

                while (idx.pos() < gap_start) {
                    idx.advance();
                }
            } while (gap_start < gap_end);
        }

        if (gap_end == n) {
            break;
        }

        factor lpf {
            .src = LPF[cur_lpf].src,
            .len = LPF[cur_lpf].end - LPF[cur_lpf].beg
        };
        
        if (idx.pos() == gap_end) {
            pos_t beg_nxt_lpf = cur_lpf + 1 == num_lpf ? n : LPF[cur_lpf + 1].beg;
            factor f = longest_prev_occ(idx, gap_end, beg_nxt_lpf - gap_end);

            if (f.len > lpf.len) {
                lpf = f;
            }
        }

        #ifndef NDEBUG
        assert(lpf.src < gap_end);

        for (pos_t i = 0; i < lpf.len; i++) {
            assert(T[gap_end + i] == T[lpf.src + i]);
        }
        #endif

        *out_it++ = lpf;
        num_phr++;
        gap_start += std::max<pos_t>(1, lpf.len);
        cur_lpf++;
    }

    if (log) {
        time = log_runtime(time);
        #ifndef NDEBUG
        std::cout << "rate of initialized values in the rolling hash index: " << idx.rate_init() << std::endl;
        #endif
    }
}