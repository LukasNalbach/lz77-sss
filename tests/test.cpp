#include <gtest/gtest.h>
#include <lz77_sss_approx/lz77_sss_approx.hpp>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<uint32_t> input_size_distrib(10000,200000);

TEST(test_lz77_sss_approx,fuzzy_test) {
    auto start_time = now();

    while (time_diff_min(start_time,now()) < 60) {
        // choose a random input length
        uint32_t input_size = input_size_distrib(gen);

        // generate a random string
        auto input = random_repetitive_string(input_size);

        // compute the factorization
        auto factorization = lz77_sss_approx<>::factorize(input);

        // decode the factorization
        auto input_decoded = lz77_sss_approx<>::decode(factorization, input.size());

        // check if the input has been decodeed correctly
        EXPECT_EQ(input,input_decoded);
    }
}