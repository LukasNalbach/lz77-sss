#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <std::input_iterator fact_it_t>
std::string lz77_sss<pos_t>::decode(fact_it_t fact_it, pos_t output_size) {
    std::string output;
    output.reserve(output_size);
    factor f;

    while (output.size() < output_size) {
        f = *fact_it++;

        if (f.len == 0) {
            output.push_back(unsigned_to_char<pos_t>(f.src));
        } else {
            #ifndef NDEBUG
            assert(f.src < output.size());
            #endif

            for (pos_t i = 0; i < f.len; i++) {
                output.push_back(output[f.src + i]);
            }
        }
    }

    return output;
}

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::build_lce() {
    if (log) {
        std::cout << "building LCE data structure" << std::flush;
    }

    uint64_t alloc_before = malloc_count_current();
    LCE = lce_t(T);
    size_sss = LCE.get_sync_set().size();
    uint64_t bytes_ssa = size_sss * sizeof(uint32_t);
    uint64_t bytes_lce = malloc_count_current() - alloc_before - bytes_ssa;

    if (log) {
        std::cout << " (" << format_size(bytes_lce) << ")";
        time = log_runtime(time);
        std::cout << "tau = " << tau;
        std::cout << ", SSA size: " << format_size(bytes_ssa) << std::endl;
        std::cout << "peak memory consumption: " << format_size(
            malloc_count_peak() - baseline_memory_alloc) << std::endl;
        std::cout << "input size / SSS size = " << n / (double) size_sss << std::endl;
    }
}