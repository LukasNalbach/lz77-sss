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
#include <lz77_sss/data_structures/rabin_karp_substring.hpp>

#include "test-progress.hpp"

thread_local std::mt19937 gen(std::random_device{}());
thread_local std::uniform_int_distribution<uint64_t> window_distrib(1, 1000);
thread_local std::uniform_int_distribution<uint64_t> roll_dist_distrib(1, 10000);
thread_local std::uniform_int_distribution<uint64_t> substr_len_distrib(1, 10000);
thread_local std::uniform_int_distribution<uint64_t> sampl_rate_distrib(1, 10000);

void test_roll()
{
    std::string input = random_repetitive_string(1, 30000);
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    uint64_t window = std::min<uint64_t>(window_distrib(gen), input.size());
    rabin_karp_substring<31, uint32_t> rks(input.data(), input.size(), sampl_rate_distrib(gen), window, num_threads);

    if (window < input.size()) {
        std::uniform_int_distribution<uint64_t> roll_pos_distrib(0, input.size() - window - 1);
        for (uint64_t i = 0; i < 100; i++) {
            uint64_t pos = roll_pos_distrib(gen);
            uint64_t dist = std::min(roll_dist_distrib(gen), input.size() - (pos + window));
            uint64_t fp = rks.substring(pos, window);
            for (uint64_t d = 0; d < dist; d++) fp = rks.roll(fp, input[pos + d], input[pos + d + window]);
            uint64_t fp_dest = rks.substring(pos + dist, window);
            EXPECT_EQ(fp, fp_dest);
        }
    }
}

void test_concat()
{
    std::string input = random_repetitive_string(1, 30000);
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    uint64_t window = std::min<uint64_t>(window_distrib(gen), input.size());
    rabin_karp_substring<31, uint32_t> rks(input.data(), input.size(), sampl_rate_distrib(gen), window, num_threads);
    std::uniform_int_distribution<uint64_t> substr_pos_distrib(0, input.size() - 1);

    for (uint64_t i = 0; i < 100; i++) {
        uint64_t pos = substr_pos_distrib(gen);
        uint64_t len = std::min(substr_len_distrib(gen), input.size() - pos);
        uint64_t len_mid = len / 2;
        uint64_t fp_left = rks.substring(pos, len_mid);
        uint64_t fp_right = rks.substring(pos + len_mid, len - len_mid);
        uint64_t concat = rks.concat(fp_left, fp_right, len - len_mid);
        uint64_t fp_full = rks.substring_naive(pos, len);
        EXPECT_EQ(fp_full, concat);
    }
}

void test_substring()
{
    std::string input = random_repetitive_string(1, 30000);
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    uint64_t window = std::min<uint64_t>(window_distrib(gen), input.size());
    rabin_karp_substring<31, uint32_t> rks(input.data(), input.size(), sampl_rate_distrib(gen), window, num_threads);
    std::uniform_int_distribution<uint64_t> substr_pos_distrib(0, input.size() - 1);

    for (uint64_t i = 0; i < 100; i++) {
        uint64_t pos = substr_pos_distrib(gen);
        uint64_t len = std::min(substr_len_distrib(gen), input.size() - pos);
        uint64_t fp_naive = rks.substring_naive(pos, len);
        uint64_t fp = rks.substring(pos, len);
        EXPECT_EQ(fp_naive, fp);
    }
}

TEST(test_rabin_karp_substring, roll)
{
    run_fuzz("rabin-karp-substring", {
        { "roll", [](uint64_t) { test_roll(); }, false },
    }, fuzz_iterations(5000));
}

TEST(test_rabin_karp_substring, concat)
{
    run_fuzz("rabin-karp-substring", {
        { "concat", [](uint64_t) { test_concat(); }, false },
    }, fuzz_iterations(5000));
}

TEST(test_rabin_karp_substring, substring)
{
    run_fuzz("rabin-karp-substring", {
        { "substring", [](uint64_t) { test_substring(); }, false },
    }, fuzz_iterations(5000));
}