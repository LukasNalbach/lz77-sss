#include <gtest/gtest.h>
#include <lz77_sss/lz77_sss.hpp>

std::random_device rd;
std::mt19937 gen(rd());

std::vector<lz77_sss<>::factor> run_alg(std::string& input) {
    switch (std::rand() % 11) {
        case 0: return lz77_sss<>::factorize_approximate<greedy_naive, lpf_all>(input);
        case 1: return lz77_sss<>::factorize_approximate<greedy, lpf_naive>(input);
        case 2: return lz77_sss<>::factorize_approximate<greedy, lpf_all>(input);
        case 3: return lz77_sss<>::factorize_approximate<greedy, lpf_lnf_all>(input);
        case 4: return lz77_sss<>::factorize_approximate<blockwise_all, lpf_all>(input);
        case 5: return lz77_sss<>::factorize_exact<greedy, lpf_all, naive, semi_dynamic_square_grid>(input);
        case 6: return lz77_sss<>::factorize_exact<greedy, lpf_all, naive, static_weighted_square_grid>(input);
        case 7: return lz77_sss<>::factorize_exact<greedy, lpf_all, with_samples, semi_dynamic_square_grid>(input);
        case 8: return lz77_sss<>::factorize_exact<greedy, lpf_all, with_samples, static_weighted_square_grid>(input);
        case 9: return lz77_sss<>::factorize_exact<greedy, lpf_all, without_samples, semi_dynamic_square_grid>(input);
        case 10: return lz77_sss<>::factorize_exact<greedy, lpf_all, without_samples, static_weighted_square_grid>(input);
    }
}

TEST(test_lz77_sss, fuzzy_test) {
    auto start_time = now();

    while (time_diff_min(start_time,now()) < 60) {
        // generate a random string
        auto input = random_repetitive_string(10000, 100000);

        // compute the factorization
        auto factorization = run_alg(input);

        // decode the factorization
        auto input_decoded = lz77_sss<>::decode(factorization, input.size());

        // check if the input has been decodeed correctly
        EXPECT_EQ(input,input_decoded);
    }
}