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
#include <lz77_sss/lz77_sss.hpp>

#include "test-progress.hpp"

std::random_device rd;
std::mt19937 gen(rd());

using factor_out = std::function<void(lz77_sss<>::factor)>;

void run_approximate(std::string& in, const factor_out& out, uint16_t num_threads) {
    switch (std::rand() % 8) {
        case 0: lz77_sss<>::factorize_approximate<greedy_naive, lpf_naive>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 1: lz77_sss<>::factorize_approximate<greedy_naive, lpf_opt>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 2: lz77_sss<>::factorize_approximate<greedy_naive, lpf_lnf_naive>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 3: lz77_sss<>::factorize_approximate<greedy_naive, lpf_lnf_opt>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 4: lz77_sss<>::factorize_approximate<greedy, lpf_naive>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 5: lz77_sss<>::factorize_approximate<greedy, lpf_opt>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 6: lz77_sss<>::factorize_approximate<greedy, lpf_lnf_naive>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 7: lz77_sss<>::factorize_approximate<greedy, lpf_lnf_opt>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
    }
}

void run_exact_naive(std::string& in, const factor_out& out, uint16_t num_threads) {
    switch (std::rand() % 4) {
        case 0: lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, semi_dynamic_square_grid>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 1: lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, static_weighted_kd_tree>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 2: lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, decomposed_semi_dynamic_square_grid>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 3: lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, decomposed_static_weighted_kd_tree>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
    }
}

void run_exact_with_samples(std::string& in, const factor_out& out, uint16_t num_threads) {
    switch (std::rand() % 2) {
        case 0: lz77_sss<>::factorize_exact<greedy, lpf_opt, with_samples, decomposed_semi_dynamic_square_grid>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 1: lz77_sss<>::factorize_exact<greedy, lpf_opt, with_samples, decomposed_static_weighted_kd_tree>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
    }
}

void run_exact_without_samples(std::string& in, const factor_out& out, uint16_t num_threads) {
    switch (std::rand() % 2) {
        case 0: lz77_sss<>::factorize_exact<greedy, lpf_opt, without_samples, decomposed_semi_dynamic_square_grid>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
        case 1: lz77_sss<>::factorize_exact<greedy, lpf_opt, without_samples, decomposed_static_weighted_kd_tree>(in.data(), in.size(), out, { .num_threads = num_threads }); break;
    }
}

void fuzz_factorize(const std::function<void(std::string&, const factor_out&, uint16_t)>& run, uint32_t max_input_size) {
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    std::string input = random_repetitive_string(1, max_input_size);
    std::vector<lz77_sss<>::factor> factorization;
    run(input, [&](lz77_sss<>::factor f) { factorization.emplace_back(f); }, num_threads);
    std::string input_decoded;
    no_init_resize(input_decoded, input.size());
    lz77_sss<>::decode(factorization.begin(), input_decoded.data(), input.size());
    EXPECT_EQ(input, input_decoded);
}

TEST(test_lz77_sss, approximate) {
    run_fuzz("lz77-sss", {
        { "approximate", [](uint64_t) { fuzz_factorize(run_approximate, 100000); }, false },
    }, fuzz_iterations(1000));
}

TEST(test_lz77_sss, exact_naive) {
    run_fuzz("lz77-sss", {
        { "exact-naive", [](uint64_t) { fuzz_factorize(run_exact_naive, 10000); }, false },
    }, fuzz_iterations(1000));
}

TEST(test_lz77_sss, exact_with_samples) {
    run_fuzz("lz77-sss", {
        { "exact-with-samples", [](uint64_t) { fuzz_factorize(run_exact_with_samples, 10000); }, false },
    }, fuzz_iterations(1000));
}

TEST(test_lz77_sss, exact_without_samples) {
    run_fuzz("lz77-sss", {
        { "exact-without-samples", [](uint64_t) { fuzz_factorize(run_exact_without_samples, 10000); }, false },
    }, fuzz_iterations(1000));
}