#include <gtest/gtest.h>
#include <lz77_sss/lz77_sss.hpp>

std::random_device rd;
std::mt19937 gen(rd());

void run_alg(std::string& input, lz77_sss<>::output_it_t&& output) {
    switch (std::rand() % 16) {
        case 0: lz77_sss<>::factorize_approximate<greedy_naive, lpf_naive>(input.data(), input.size(), output);
        case 1: lz77_sss<>::factorize_approximate<greedy_naive, lpf_opt>(input.data(), input.size(), output);
        case 2: lz77_sss<>::factorize_approximate<greedy_naive, lpf_lnf_naive>(input.data(), input.size(), output);
        case 3: lz77_sss<>::factorize_approximate<greedy_naive, lpf_lnf_opt>(input.data(), input.size(), output);
        case 4: lz77_sss<>::factorize_approximate<greedy, lpf_naive>(input.data(), input.size(), output);
        case 5: lz77_sss<>::factorize_approximate<greedy, lpf_opt>(input.data(), input.size(), output);
        case 6: lz77_sss<>::factorize_approximate<greedy, lpf_lnf_naive>(input.data(), input.size(), output);
        case 7: lz77_sss<>::factorize_approximate<greedy, lpf_lnf_opt>(input.data(), input.size(), output);
        case 8: lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, semi_dynamic_square_grid>(input.data(), input.size(), output);
        case 9: lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, static_weighted_kd_tree>(input.data(), input.size(), output);
        case 10: lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, decomposed_semi_dynamic_square_grid>(input.data(), input.size(), output);
        case 11: lz77_sss<>::factorize_exact<greedy, lpf_opt, naive, decomposed_static_weighted_kd_tree>(input.data(), input.size(), output);
        case 12: lz77_sss<>::factorize_exact<greedy, lpf_opt, with_samples, decomposed_semi_dynamic_square_grid>(input.data(), input.size(), output);
        case 13: lz77_sss<>::factorize_exact<greedy, lpf_opt, with_samples, decomposed_static_weighted_kd_tree>(input.data(), input.size(), output);
        case 14: lz77_sss<>::factorize_exact<greedy, lpf_opt, without_samples, decomposed_semi_dynamic_square_grid>(input.data(), input.size(), output);
        case 15: lz77_sss<>::factorize_exact<greedy, lpf_opt, without_samples, decomposed_static_weighted_kd_tree>(input.data(), input.size(), output);
    }
}

TEST(test_lz77_sss, fuzzy_test) {
    auto start_time = now();

    while (time_diff_min(start_time,now()) < 60) {
        // generate a random string
        std::string input = random_repetitive_string(4'000, 1'000'000);

        // compute the factorization
        std::vector<lz77_sss<>::factor> factorization;
        run_alg(input, [&](auto f){factorization.emplace_back(f);});

        // decode the factorization
        std::string input_decoded;
        no_init_resize(input_decoded, input.size());
        lz77_sss<>::decode(factorization.begin(), input_decoded.data(), input.size());

        // check if the input has been decodeed correctly
        EXPECT_EQ(input, input_decoded);
    }
}