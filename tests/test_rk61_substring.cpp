#include <gtest/gtest.h>
#include <lz77_sss/data_structures/rk61_substring.hpp>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<uint64_t> window_distrib(1, 1000);
std::uniform_int_distribution<uint64_t> roll_dist_distrib(1, 10000);
std::uniform_int_distribution<uint64_t> substr_len_distrib(1, 10000);
std::uniform_int_distribution<uint64_t> sampl_rate_distrib(1, 10000);

std::string input;

TEST(test_rk_substr, fuzzy_test)
{
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        uint64_t window = window_distrib(gen);

        // choose a random input
        input = random_repetitive_string(10000, 200000);
        std::uniform_int_distribution<uint64_t> substr_pos_distrib(
            0, input.size() - window - 1);

        // build the data structure
        rk61_substring rks(input, sampl_rate_distrib(gen), window, omp_get_max_threads());

        // check if rolling fingerprints work
        for (uint64_t i = 0; i < 1000; i++) {
            uint64_t pos = substr_pos_distrib(gen);
            uint64_t dist = std::min(roll_dist_distrib(gen), input.size() - (pos + window));
            uint64_t fp = rks.substring<>(pos, window);
            for (uint64_t d = 0; d < dist; d++) fp = rks.roll(fp, input[pos + d], input[pos + d + window]);
            uint64_t fp_dest = rks.substring<>(pos + dist, window);
            EXPECT_EQ(fp, fp_dest);
        }

        // check if fingerprints are concatenated correctly
        for (uint64_t i = 0; i < 1000; i++) {
            uint64_t pos = substr_pos_distrib(gen);
            uint64_t len = std::min(substr_len_distrib(gen), input.size() - pos);
            uint64_t len_mid = len / 2;
            uint64_t fp_left = rks.substring<>(pos, len_mid);
            uint64_t fp_right = rks.substring<>(pos + len_mid, len - len_mid);
            uint64_t concat = rks.concat(fp_left, fp_right, len - len_mid);
            uint64_t fp_full = rks.substring_naive<>(pos, len);
            EXPECT_EQ(fp_full, concat);
        }

        // extract fingerprints of random substrings and check if they are correct
        for (uint64_t i = 0; i < 1000; i++) {
            uint64_t pos = substr_pos_distrib(gen);
            uint64_t len = std::min(substr_len_distrib(gen), input.size() - pos);
            uint64_t fp_naive = rks.substring_naive<>(pos, len);
            uint64_t fp = rks.substring<>(pos, len);
            EXPECT_EQ(fp_naive, fp);
        }
    }
}