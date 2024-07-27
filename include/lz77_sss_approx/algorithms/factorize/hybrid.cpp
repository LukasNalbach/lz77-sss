#pragma once

#include <lz77_sss_approx/lz77_sss_approx.hpp>
#include <list>

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
template <pos_t... patt_lens>
lz77_sss_approx<pos_t>::lpf
lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::lpf_extend_left(
    rolling_hash_index_107<pos_t, patt_lens...>& idx, pos_t pos, pos_t max_lce_l
) {
    assert(idx.pos() == pos);
    lpf phr{0, 0, 0};

    for_constexpr<sizeof...(patt_lens) - 1, - 1, - 1>([&](auto patt_len_idx) {
        if (phr.end == phr.beg) {
            pos_t src = idx.template advance_and_get_occ<patt_len_idx>();

            if (src < pos) {
                pos_t lce_l = src == 0 ? 0 : LCE_L(src - 1, pos - 1, max_lce_l);
                pos_t end = pos + LCE_R(src, pos);
                pos_t beg = pos - lce_l;

                if (end - beg > phr.end - phr.beg) {
                    phr.beg = beg;
                    phr.end = end;
                    phr.src = src - lce_l;

                    #ifndef NDEBUG
                    assert(src < pos);

                    for (pos_t j = 0; j < end - beg; j++) {
                        assert(T[phr.src + j] == T[phr.beg + j]);
                    }
                    #endif
                }
            }
        } else {
            idx.template advance<patt_len_idx>();
        }
    });

    idx.inc_pos();

    return phr;
}

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
template <pos_t... patt_lens>
void lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::factorize_hybrid() {
    if (log) {
        std::cout << "initializing compact index" << std::flush;
    }

    rolling_hash_index_107<pos_t, patt_lens...> idx(T, target_index_size);

    if (log) {
        std::cout << " (size: " << format_size(idx.size_in_bytes()) << ")";
        time = log_runtime(time);
        std::cout << "factorizing" << std::flush;
    }

    uint64_t max_num_phr = 6;

    std::list<lpf> P;
    typename std::list<lpf>::iterator phr_m = P.end();
    typename std::list<lpf>::reverse_iterator phr_l = P.rend();
    typename std::list<lpf>::iterator phr_r = P.end();

    pos_t x = 0;
    pos_t i = 0;
    pos_t j = 0;

    while (i < n) {
        lpf phr = lpf_extend_left(idx, i, i - j);

        if (phr.end - phr.beg != 0) {
            if (P.empty()) {
                P.emplace_back(phr);
                phr_m = P.begin();
            } else {
                phr_l = --std::make_reverse_iterator(phr_m);
                while (phr_l != P.rend() && (*phr_l).beg >= phr.beg - 1) phr_l++;

                phr_r = phr_m;
                while (phr_r != P.end() && (*phr_r).end <= phr.end + 1) phr_r++;

                if (
                    phr_l == P.rend() || phr_r == P.end() ||
                    (*phr_l).beg != (*phr_r).beg
                ) {
                    if (phr_l != P.rend() && (*phr_l).end > phr.beg) {
                        (*phr_l).end = phr.beg;
                    }

                    if (phr_r != P.end() && (*phr_r).beg < phr.end) {
                        (*phr_r).src += phr.end - (*phr_r).beg;
                        (*phr_r).beg = phr.end;
                    }

                    phr_m = P.emplace(
                        P.erase(phr_l.base(), phr_r),
                        phr
                    );

                    while (P.size() > max_num_phr) {
                        while (j < P.front().beg) {
                            *out_it++ = factor{char_to_uchar(T[j]), 0};
                            num_phr++;
                            j++;
                        }

                        *out_it++ = factor{
                            .src = P.front().src,
                            .len = P.front().end - P.front().beg
                        };

                        #ifndef NDEBUG
                        assert(j == P.front().beg);
                        assert(P.front().src < j);

                        for (pos_t j = 0; j < P.front().end - P.front().beg; j++) {
                            assert(T[P.front().src + j] == T[P.front().beg + j]);
                        }
                        #endif

                        j = P.front().end;
                        P.pop_front();
                    }
                }
            }
        }

        i = std::max<pos_t>(i + 1, phr.end);

        while (idx.pos() < i) {
            idx.inc_pos();
        }
    }

    while (!P.empty()) {
        while (j < P.front().beg) {
            *out_it++ = factor{char_to_uchar(T[j]), 0};
            num_phr++;
            j++;
        }

        *out_it++ = factor{
            .src = P.front().src,
            .len = P.front().end - P.front().beg
        };

        #ifndef NDEBUG
        assert(j == P.front().beg);
        assert(P.front().src < j);

        for (pos_t j = 0; j < P.front().end - P.front().beg; j++) {
            assert(T[P.front().src + j] == T[P.front().beg + j]);
        }
        #endif

        j = P.front().end;
        P.pop_front();
    }

    while (j < n) {
        *out_it++ = factor{char_to_uchar(T[j]), 0};
        num_phr++;
        j++;
    }

    assert(idx.pos() == n);

    if (log) {
        time = log_runtime(time);
        #ifndef NDEBUG
        std::cout << "rate of initialized values in the compact index: " << idx.rate_init() << std::endl;
        #endif
    }
}