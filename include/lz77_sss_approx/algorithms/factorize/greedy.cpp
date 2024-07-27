#pragma once

#include <lz77_sss_approx/lz77_sss_approx.hpp>

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
template <pos_t... patt_lens>
lz77_sss_approx<pos_t>::factor
lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::longest_prev_occ(
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
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
template <pos_t... patt_lens>
void lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::factorize_greedy() {
    if (log) {
        std::cout << "initializing compact index" << std::flush;
    }

    rolling_hash_index_107<pos_t, patt_lens...> idx(T, target_index_size);

    if (log) {
        std::cout << " (size: " << format_size(idx.size_in_bytes()) << ")";
        time = log_runtime(time);
        std::cout << "factorizing" << std::flush;
    }
    
    pos_t cur_lpf = 0;
    pos_t gap_start = 1;

    if (num_lpf != 0 && LPF[0].beg == 0) {
        cur_lpf = 1;
    }

    *out_it++ = factor{char_to_uchar(T[0]),0};
    num_phr++;
    idx.advance();
    
    while (true) {
        pos_t gap_end = cur_lpf == num_lpf ? n : LPF[cur_lpf].beg;

        while (gap_start < gap_end) {
            factor f = longest_prev_occ(idx, gap_start, gap_end - gap_start);
            *out_it++ = f;
            num_phr++;
            gap_start += std::max<pos_t>(1, f.len);

            while (idx.pos() < gap_start) {
                idx.advance();
            }
        }

        if (gap_end == n) {
            break;
        }

        factor lpf {
            .src = LPF[cur_lpf].src,
            .len = LPF[cur_lpf].end - LPF[cur_lpf].beg
        };

        #ifndef NDEBUG
        assert(lpf.src < gap_end);

        for (pos_t i = 0; i < lpf.len; i++) {
            assert(T[gap_end + i] == T[lpf.src + i]);
        }
        #endif

        {
            pos_t beg_nxt_lpf = cur_lpf + 1 == num_lpf ? n : LPF[cur_lpf + 1].beg;
            factor f = longest_prev_occ(idx, gap_end, beg_nxt_lpf - gap_end);

            if (f.len > lpf.len) {
                lpf = f;
            }
        }

        *out_it++ = lpf;
        num_phr++;
        gap_start += std::max<pos_t>(1, lpf.len);

        idx.reinit(gap_start);
        /*
        while (idx.pos() < gap_start) {
            idx.advance();
        }
        */

        cur_lpf++;
    }

    if (log) {
        time = log_runtime(time);
        #ifndef NDEBUG
        std::cout << "rate of initialized values in the compact index: " << idx.rate_init() << std::endl;
        #endif
    }
}