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
#include <text/packed_text.hpp>
#include <text/split_text.hpp>
#include <lz77_sss/misc/repetitive_input.hpp>

#include "test-progress.hpp"

thread_local std::mt19937 gen(std::random_device{}());
thread_local std::uniform_int_distribution<uint64_t> window_distrib(1, 1000);
thread_local std::uniform_int_distribution<uint64_t> roll_dist_distrib(1, 10000);
thread_local std::uniform_int_distribution<uint64_t> substr_len_distrib(1, 10000);
thread_local std::uniform_int_distribution<uint64_t> sample_rate_distrib(1, 10000);

template <uint8_t mers_exp>
void test_roll()
{
    std::string input = random_repetitive_input<std::string>(1, 30000);
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    uint64_t window = std::min<uint64_t>(window_distrib(gen), input.size());
    rabin_karp_substring<mers_exp> rks(input.data(), input.size(), sample_rate_distrib(gen), window, num_threads);

    if (window < input.size()) {
        std::uniform_int_distribution<uint64_t> roll_pos_distrib(0, input.size() - window - 1);
        for (uint64_t i = 0; i < 100; i++) {
            uint64_t pos = roll_pos_distrib(gen);
            uint64_t dist = std::min(roll_dist_distrib(gen), input.size() - (pos + window));
            uint64_t fp = rks.substr_fp(pos, window);
            for (uint64_t d = 0; d < dist; d++) fp = rks.roll(fp, input[pos + d], input[pos + d + window]);
            uint64_t fp_dest = rks.substr_fp(pos + dist, window);
            EXPECT_EQ(fp, fp_dest);
        }
    }
}

template <uint8_t mers_exp>
void test_concat()
{
    std::string input = random_repetitive_input<std::string>(1, 30000);
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    uint64_t window = std::min<uint64_t>(window_distrib(gen), input.size());
    rabin_karp_substring<mers_exp> rks(input.data(), input.size(), sample_rate_distrib(gen), window, num_threads);
    std::uniform_int_distribution<uint64_t> substr_pos_distrib(0, input.size() - 1);

    for (uint64_t i = 0; i < 100; i++) {
        uint64_t pos = substr_pos_distrib(gen);
        uint64_t len = std::min(substr_len_distrib(gen), input.size() - pos);
        uint64_t len_mid = len / 2;
        uint64_t fp_left = rks.substr_fp(pos, len_mid);
        uint64_t fp_right = rks.substr_fp(pos + len_mid, len - len_mid);
        uint64_t concat = rks.concat(fp_left, fp_right, len - len_mid);
        uint64_t fp_full = rks.substr_fp_naive(pos, len);
        EXPECT_EQ(fp_full, concat);
    }
}

template <uint8_t mers_exp>
void test_substring()
{
    std::string input = random_repetitive_input<std::string>(1, 30000);
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    uint64_t window = std::min<uint64_t>(window_distrib(gen), input.size());
    rabin_karp_substring<mers_exp> rks(input.data(), input.size(), sample_rate_distrib(gen), window, num_threads);
    std::uniform_int_distribution<uint64_t> substr_pos_distrib(0, input.size() - 1);

    for (uint64_t i = 0; i < 100; i++) {
        uint64_t pos = substr_pos_distrib(gen);
        uint64_t len = std::min(substr_len_distrib(gen), input.size() - pos);
        uint64_t fp_naive = rks.substr_fp_naive(pos, len);
        uint64_t fp = rks.substr_fp(pos, len);
        EXPECT_EQ(fp_naive, fp);
    }
}

template <typename fnc_t>
void with_random_int_text(fnc_t fnc)
{
    const uint32_t sigma = random_log_uniform_size(1, uint64_t { 1 } << 16, gen);
    std::vector<uint32_t> input = random_repetitive_input<std::vector<uint32_t>>(1, 10000, uint32_t(0), uint32_t(sigma - 1));
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);

    switch (std::uniform_int_distribution<int>(0, 2)(gen)) {
        case 0: fnc(lce::text::direct_text<uint32_t>(input.data(), input.size(), sigma), num_threads); break;
        case 1: fnc(lce::text::packed_text<uint32_t>(input.data(), input.size(), sigma, num_threads), num_threads); break;
        default: fnc(lce::text::split_text<uint32_t>(input.data(), input.size(), sigma, num_threads), num_threads); break;
    }
}

void test_roll_int()
{
    with_random_int_text([](const auto& text, uint16_t num_threads) {
        using text_t = std::decay_t<decltype(text)>;
        const uint64_t n = text.size();
        uint64_t window = std::min<uint64_t>(window_distrib(gen), n);
        rabin_karp_substring<61, text_t> rks(text, n, sample_rate_distrib(gen), window, num_threads);

        if (window < n) {
            std::uniform_int_distribution<uint64_t> roll_pos_distrib(0, n - window - 1);
            for (uint64_t i = 0; i < 100; i++) {
                uint64_t pos = roll_pos_distrib(gen);
                uint64_t dist = std::min(roll_dist_distrib(gen) / 10, n - (pos + window));
                uint64_t fp = rks.substr_fp(pos, window);
                for (uint64_t d = 0; d < dist; d++) fp = rks.roll(fp, text[pos + d], text[pos + d + window]);
                uint64_t fp_dest = rks.substr_fp(pos + dist, window);
                EXPECT_EQ(fp, fp_dest);
            }
        }
    });
}

void test_concat_int()
{
    with_random_int_text([](const auto& text, uint16_t num_threads) {
        using text_t = std::decay_t<decltype(text)>;
        const uint64_t n = text.size();
        uint64_t window = std::min<uint64_t>(window_distrib(gen), n);
        rabin_karp_substring<61, text_t> rks(text, n, sample_rate_distrib(gen), window, num_threads);
        std::uniform_int_distribution<uint64_t> substr_pos_distrib(0, n - 1);

        for (uint64_t i = 0; i < 100; i++) {
            uint64_t pos = substr_pos_distrib(gen);
            uint64_t len = std::min(substr_len_distrib(gen), n - pos);
            uint64_t len_mid = len / 2;
            uint64_t fp_left = rks.substr_fp(pos, len_mid);
            uint64_t fp_right = rks.substr_fp(pos + len_mid, len - len_mid);
            uint64_t concat = rks.concat(fp_left, fp_right, len - len_mid);
            uint64_t fp_full = rks.substr_fp_naive(pos, len);
            EXPECT_EQ(fp_full, concat);
        }
    });
}

void test_substring_int()
{
    with_random_int_text([](const auto& text, uint16_t num_threads) {
        using text_t = std::decay_t<decltype(text)>;
        const uint64_t n = text.size();
        uint64_t window = std::min<uint64_t>(window_distrib(gen), n);
        rabin_karp_substring<61, text_t> rks(text, n, sample_rate_distrib(gen), window, num_threads);
        std::uniform_int_distribution<uint64_t> substr_pos_distrib(0, n - 1);

        for (uint64_t i = 0; i < 100; i++) {
            uint64_t pos = substr_pos_distrib(gen);
            uint64_t len = std::min(substr_len_distrib(gen), n - pos);
            uint64_t fp_naive = rks.substr_fp_naive(pos, len);
            uint64_t fp = rks.substr_fp(pos, len);
            EXPECT_EQ(fp_naive, fp);

            uint64_t other = substr_pos_distrib(gen);

            if (other + len <= n) {
                bool same = true;
                for (uint64_t k = 0; k < len && same; k++) same = text[pos + k] == text[other + k];
                EXPECT_EQ(same, fp == rks.substr_fp(other, len));
            }
        }
    });
}

TEST(test_rabin_karp_substring, roll_mersenne_31)
{
    run_fuzz("rabin-karp-substring", {
        { "roll-mersenne-31", [](uint64_t) { test_roll<31>(); }, false },
    }, fuzz_iterations(2500));
}

TEST(test_rabin_karp_substring, roll_mersenne_61)
{
    run_fuzz("rabin-karp-substring", {
        { "roll-mersenne-61", [](uint64_t) { test_roll<61>(); }, false },
    }, fuzz_iterations(2500));
}

TEST(test_rabin_karp_substring, concat_mersenne_31)
{
    run_fuzz("rabin-karp-substring", {
        { "concat-mersenne-31", [](uint64_t) { test_concat<31>(); }, false },
    }, fuzz_iterations(2500));
}

TEST(test_rabin_karp_substring, concat_mersenne_61)
{
    run_fuzz("rabin-karp-substring", {
        { "concat-mersenne-61", [](uint64_t) { test_concat<61>(); }, false },
    }, fuzz_iterations(2500));
}

TEST(test_rabin_karp_substring, substring_mersenne_31)
{
    run_fuzz("rabin-karp-substring", {
        { "substring-mersenne-31", [](uint64_t) { test_substring<31>(); }, false },
    }, fuzz_iterations(2500));
}

TEST(test_rabin_karp_substring, substring_mersenne_61)
{
    run_fuzz("rabin-karp-substring", {
        { "substring-mersenne-61", [](uint64_t) { test_substring<61>(); }, false },
    }, fuzz_iterations(2500));
}

TEST(test_rabin_karp_substring, roll_int)
{
    run_fuzz("rabin-karp-substring", {
        { "roll-int", [](uint64_t) { test_roll_int(); }, false },
    }, fuzz_iterations(1500));
}

TEST(test_rabin_karp_substring, concat_int)
{
    run_fuzz("rabin-karp-substring", {
        { "concat-int", [](uint64_t) { test_concat_int(); }, false },
    }, fuzz_iterations(1500));
}

TEST(test_rabin_karp_substring, substring_int)
{
    run_fuzz("rabin-karp-substring", {
        { "substring-int", [](uint64_t) { test_substring_int(); }, false },
    }, fuzz_iterations(1500));
}
