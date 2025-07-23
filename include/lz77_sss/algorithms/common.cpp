#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <std::input_iterator fact_it_t, typename char_t>
void lz77_sss<pos_t>::decode(fact_it_t fact_it, char_t* output, pos_t output_size)
{
    factor f;
    char_t* output_it = output;
    char_t* end = output + output_size;

    while (output_it < end) {
        f = *fact_it++;

        if (f.len == 0) {
            *output_it++ = unsigned_to_char<pos_t, char_t>(f.src);
        } else {
            #ifndef NDEBUG
            assert(f.src < output_size);
            #endif

            for (pos_t i = 0; i < f.len; i++) {
                *output_it++ = output[f.src + i];
            }
        }
    }
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::build_lce()
{
    if (log) {
        std::cout << "building LCE data structure" << std::flush;
    }

    uint64_t alloc_before = malloc_count_current();
    LCE = lce_t(T, n);
    size_sss = LCE.get_sync_set().size();
    uint64_t bytes_sa_s = size_sss * sizeof(uint32_t);
    uint64_t bytes_lce = malloc_count_current() - alloc_before - bytes_sa_s;

    if (log) {
        std::string phase = LPF[0].empty() ? "lce_t" : "lce_rev_t";
        log_phase(phase, time_diff_ns(time, now()));
        std::cout << " (" << format_size(bytes_lce) << ")";
        time = log_runtime(time);
        std::cout << "tau = " << tau;
        std::cout << ", SA_S size: " << format_size(bytes_sa_s) << std::endl;
        std::cout << "peak memory consumption: " << format_size(malloc_count_peak() - baseline_memory_alloc) << std::endl;
        std::cout << "input length / SSS size = " << n / (double)size_sss << std::endl;
    }
}