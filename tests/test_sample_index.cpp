/**
 * part of LukasNalbach/lz77-sss
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <gtest/gtest.h>
#include <ips4o.hpp>
#include <lz77_sss/data_structures/sample_index/sample_index.hpp>
#include <text/packed_text.hpp>
#include <text/split_text.hpp>
#include <lz77_sss/misc/repetitive_input.hpp>

#include "test-progress.hpp"

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<uint32_t> avg_sample_rate_distrib(1, 10);
std::uniform_int_distribution<uint32_t> pattern_length_distrib(1, 1000);
std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);

template <direction dir, typename input_t, typename index_t>
void test_query(const input_t& input, const std::vector<uint40_t>& sampling, index_t& index)
{
    std::vector<uint64_t> occurrences;
    std::vector<uint64_t> correct_occurrences;
    std::uniform_int_distribution<uint32_t> pattern_pos_distrib(0, input.size() - 1);
    uint32_t pattern_pos = pattern_pos_distrib(gen);
    typename index_t::query_ctx_t query = index.query();
    uint32_t pattern_length = std::min<uint32_t>(
        dir == LEFT ? pattern_pos + 1 : (input.size() - pattern_pos),
        pattern_length_distrib(gen));
    std::uniform_int_distribution<uint32_t> step_size_distrib(1, 3);
    uint32_t current_pattern_length = 0;
    bool occurs = true;

    while (occurs && current_pattern_length < pattern_length) {
        current_pattern_length = std::min<uint32_t>(
            pattern_length, current_pattern_length + step_size_distrib(gen));
        occurs = index.template extend<dir>(query, pattern_pos, current_pattern_length);
    }

    if (occurs) {
        index.template locate<dir>(query, occurrences);
        ips4o::sort(occurrences.begin(), occurrences.end());
    }

    for (uint32_t i = 0; i < sampling.size(); i++) {
        uint32_t sample_pos = sampling[i];
        occurs = true;

        if (dir == LEFT ?
            (sample_pos >= pattern_length - 1) :
            (sample_pos + pattern_length <= input.size())
        ) {
            for (uint32_t j = 0; j < pattern_length; j++) {
                if (dir == LEFT ?
                    input[sample_pos - j] != input[pattern_pos - j] :
                    input[sample_pos + j] != input[pattern_pos + j]
                ) {
                    occurs = false;
                    break;
                }
            }
            if (occurs) {
                correct_occurrences.emplace_back(sample_pos);
            }
        }
    }

    EXPECT_EQ(occurrences, correct_occurrences);
}

template <typename text_t, typename input_t>
void test_extend_locate_on(const text_t& text, input_t& input, interval_samples mode, uint16_t num_threads)
{
    using sym_t = typename input_t::value_type;
    using lce_r_t = lce::ds::lce_naive_wordwise_xor<sym_t>;

    // choose a random average sample rate
    uint32_t avg_sample_rate = avg_sample_rate_distrib(gen);
    std::uniform_int_distribution<uint32_t> sample_distance_distrib(1, 2 * avg_sample_rate);

    // compute a random sampling of text positions
    std::vector<uint40_t> sampling;
    sampling.emplace_back(std::min<uint64_t>(
        input.size() - 1, sample_distance_distrib(gen)));

    while (uint64_t(sampling.back()) + 2 * avg_sample_rate < input.size()) {
        sampling.emplace_back(uint64_t(sampling.back()) + sample_distance_distrib(gen));
    }

    // build the sample-index
    sample_index<text_t, lce_r_t> index;
    const lce_r_t lce_r(input);
    index.build(text, input.size(), sampling, lce_r, mode, 64, num_threads);

    // perform random queries and check their correctness
    for (uint32_t i = 0; i < 1000; i++) {
        if (prob_distrib(gen) < 0.5) {
            test_query<LEFT>(input, sampling, index);
        } else {
            test_query<RIGHT>(input, sampling, index);
        }
    }
}

void test_extend_locate(interval_samples mode)
{
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);

    // generate a random string
    std::string input = random_repetitive_input<std::string>(1, 10000);
    test_extend_locate_on(lce::text::direct_text<char>(input.data(), input.size()), input, mode, num_threads);
}

void test_extend_locate_int(interval_samples mode)
{
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    const uint32_t sigma = random_log_uniform_size(1, uint64_t { 1 } << 18, gen);
    std::vector<uint32_t> input = random_repetitive_input<std::vector<uint32_t>>(1, 10000, uint32_t(0), uint32_t(sigma - 1));

    switch (std::uniform_int_distribution<int>(0, 2)(gen)) {
        case 0: test_extend_locate_on(lce::text::direct_text<uint32_t>(input.data(), input.size(), sigma), input, mode, num_threads); break;
        case 1: test_extend_locate_on(lce::text::packed_text<uint32_t>(input.data(), input.size(), sigma, num_threads), input, mode, num_threads); break;
        default: test_extend_locate_on(lce::text::split_text<uint32_t>(input.data(), input.size(), sigma, num_threads), input, mode, num_threads); break;
    }
}

TEST(test_sample_index, with_interval_samples)
{
    run_fuzz("sample-index", {
        { "with-interval-samples", [](uint64_t) { test_extend_locate(interval_samples::use); }, false },
    }, fuzz_iterations(1000));
}

TEST(test_sample_index, without_interval_samples)
{
    run_fuzz("sample-index", {
        { "without-interval-samples", [](uint64_t) { test_extend_locate(interval_samples::skip); }, false },
    }, fuzz_iterations(1400));
}

TEST(test_sample_index, with_interval_samples_int)
{
    run_fuzz("sample-index", {
        { "with-interval-samples-int", [](uint64_t) { test_extend_locate_int(interval_samples::use); }, false },
    }, fuzz_iterations(1000));
}

TEST(test_sample_index, without_interval_samples_int)
{
    run_fuzz("sample-index", {
        { "without-interval-samples-int", [](uint64_t) { test_extend_locate_int(interval_samples::skip); }, false },
    }, fuzz_iterations(1400));
}
