#include <gtest/gtest.h>
#include <lz77_sss/lz77_sss.hpp>

std::random_device rd;
std::mt19937 gen(rd());

std::vector<lz77_sss<>::factor> run_alg(std::string& input) {
    switch (std::rand() % 16) {
        case 0: return lz77_sss<>::factorize_approximate<greedy_naive, lpf_naive>(input);
        case 1: return lz77_sss<>::factorize_approximate<greedy_naive, lpf_opt>(input);
        case 2: return lz77_sss<>::factorize_approximate<greedy_naive, lpf_lnf_naive>(input);
        case 3: return lz77_sss<>::factorize_approximate<greedy_naive, lpf_lnf_opt>(input);
        case 4: return lz77_sss<>::factorize_approximate<greedy, lpf_naive>(input);
        case 5: return lz77_sss<>::factorize_approximate<greedy, lpf_opt>(input);
        case 6: return lz77_sss<>::factorize_approximate<greedy, lpf_lnf_naive>(input);
        case 7: return lz77_sss<>::factorize_approximate<greedy, lpf_lnf_opt>(input);
        case 8: return lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, semi_dynamic_square_grid>(input);
        case 9: return lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, static_weighted_kd_tree>(input);
        case 10: return lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, decomposed_semi_dynamic_square_grid>(input);
        case 11: return lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, decomposed_static_weighted_kd_tree>(input);
        case 12: return lz77_sss<>::factorize_exact<greedy, lpf_opt, with_samples, decomposed_semi_dynamic_square_grid>(input);
        case 13: return lz77_sss<>::factorize_exact<greedy, lpf_opt, with_samples, decomposed_static_weighted_kd_tree>(input);
        case 14: return lz77_sss<>::factorize_exact<greedy, lpf_opt, without_samples, decomposed_semi_dynamic_square_grid>(input);
        case 15: return lz77_sss<>::factorize_exact<greedy, lpf_opt, without_samples, decomposed_static_weighted_kd_tree>(input);
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