#pragma once

#include <lz77_sss_approx/lz77_sss_approx.hpp>

template <typename pos_t>
template <std::input_iterator fact_it_t>
std::string lz77_sss_approx<pos_t>::decode(fact_it_t fact_it, pos_t output_size) {
    std::string output;
    output.reserve(output_size);
    factor f;

    while (output.size() < output_size) {
        f = *fact_it++;

        if (f.len == 0) {
            output.push_back(unsigned_to_char<pos_t>(f.src));
        } else {
            for (pos_t i = 0; i < f.len; i++) {
                output.push_back(output[f.src + i]);
            }
        }
    }

    return output;
}

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
void lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::compute_char_freq() {
    if (log) {
        std::cout << "computing character frequencies in T" << std::flush;
    }

    char_freq.resize(256,0);

    for (pos_t i = 0; i < n; i++) {
        char_freq[char_to_uchar(T[i])]++;
    }

    for (uint16_t c = 0; c < 256; c++) {
        if (char_freq[c] > 0) {
            sigma++;
        }
    }

    if (log) {
        time = log_runtime(time);
        std::cout << "sigma: " << std::to_string(sigma) << std::endl;
    }
}

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
void lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::map_to_effective_alphabet() {
    if (log) {
        std::cout << "mapping T to its effective alphabet" << std::flush;
    }

    word_width = std::ceil(std::log2(sigma));
    std::vector<uint8_t> map_int(256, 0);
    map_ext.resize(256, 0);
    uint8_t c = 0;

    for (uint16_t i = 0; i < 256; i++) {
        if (char_freq[i] > 0) {
            map_int[i] = c;
            map_ext[c] = i;
            c++;
        }
    }

    for (pos_t i = 0; i < n; i++) {
        T[i] = uchar_to_char(map_int[char_to_uchar(T[i])]);
    }

    if (log) {
        time = log_runtime(time);
        std::cout << "word width: " << std::to_string(word_width) << std::endl;
    }
}

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
void lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::remap_to_original_alphabet() {
    if (log) {
        std::cout << "remapping T to its original alphabet" << std::flush;
    }

    for (pos_t i = 0; i < n; i++) {
        T[i] = uchar_to_char(map_ext[char_to_uchar(T[i])]);
    }

    if (log) {
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
void lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::build_lce() {
    if (log) {
        std::cout << "building LCE data structure" << std::flush;
    }

    LCE = lce_t(T);
    size_sss = LCE.get_sync_set().size();

    if (log) {
        time = log_runtime(time);
        std::cout << "tau = " << tau << std::endl;
        std::cout << "input size / SSS size = " << n / (double) size_sss << std::endl;
    }
}

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
pos_t lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::LCE_L(pos_t i, pos_t j, pos_t max_lce) {
    if (T[i] != T[j]) {
        return 0;
    }

    max_lce = std::min<pos_t>({i, j, max_lce});

    if (max_lce < 16) {
        pos_t min__ = i - max_lce;
        pos_t i__ = i;
        pos_t j__ = j;

        while (i__ > min__ && T[i__ - 1] == T[j__ - 1]) {
            i__--;
            j__--;
        }

        return i - i__;
    }

    uint128_t* ip1_ = reinterpret_cast<uint128_t*>(&T[i + 1]);
    uint128_t* i_ = ip1_ - 1;
    uint128_t* j_ = reinterpret_cast<uint128_t*>(&T[j + 1]) - 1;

    pos_t lce = std::countl_zero(*i_ ^ *j_) / 8;

    if (lce < 16) {
        return lce;
    }

    uint128_t* min_ = ip1_ - (max_lce / 16);

    while (i_ > min_ && *(i_ - 1) == *(j_ - 1)) {
        i_--;
        j_--;
    }

    if (i_ > min_) {
        return 16 * (ip1_ - i_) + std::countl_zero(*(i_ - 1) ^ *(j_ - 1)) / 8;
    } else {
        pos_t min__ = i - max_lce;
        pos_t i__ = reinterpret_cast<char*>(i_) - &T[0];
        pos_t j__ = reinterpret_cast<char*>(j_) - &T[0];

        while (i__ > min__ && T[i__ - 1] == T[j__ - 1]) {
            i__--;
            j__--;
        }

        return i - i__;
    }
}

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
void lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::greedy_phrase_selection(std::vector<lpf>& P) {
    ips4o::sort(P.begin(), P.end(), [](auto phr_i, auto phr_j) {
        return phr_i.beg < phr_j.beg ||
              (phr_i.beg == phr_j.beg && phr_i.end > phr_j.end);
    });

    pos_t k = 0;
    pos_t i = 1;
    pos_t p = P.size();

    while (i < p && P[i].end < P[k].end) i++;

    while (i < p) {
        pos_t x = p;

        if (i + 1 < p) {
            x = i + 1;
            
            while (x < p && P[x].beg <= P[k].end) {
                if (P[x].end > P[i].end) {
                    i = x;
                }

                x++;
            }

            if (P[i].end <= P[k].end) {
                i = x;
            }
        }

        if (i == p) {
            break;
        }

        if (P[i].beg < P[k].end) {
            P[k].end = P[i].beg;
        }

        if (P[k].end > P[k].beg) {
            k++;
        }
        
        P[k] = P[i];
        i = x;
    }
    
    P.resize(k+1);

    #ifndef NDEBUG
    for (uint32_t i = 1; i < P.size(); i++) {
        assert(P[i - 1].end <= P[i].beg);
        assert(P[i - 1].beg < P[i].beg);
    }
    #endif
}

template <typename pos_t>
template <factozize_mode fact_mode, phrase_mode phr_mode, uint64_t tau, typename out_it_t>
void lz77_sss_approx<pos_t>::factorizer<fact_mode, phr_mode, tau, out_it_t>::get_phrase_info() {
    num_lpf = LPF.size();

    if (num_lpf > 0) {
        len_lpf_phr = LPF[0].end - LPF[0].beg;
        if (LPF[0].beg > 0) num_gaps = 1;

        for (uint32_t i = 1; i < num_lpf; i++) {
            len_lpf_phr += LPF[i].end - LPF[i].beg;
            if (LPF[i].beg > LPF[i - 1].end) num_gaps++;
        }

        if (LPF[num_lpf - 1].end < n) num_gaps++;
    } else {
        num_gaps = 1;
    }

    len_gaps = n - len_lpf_phr;
}