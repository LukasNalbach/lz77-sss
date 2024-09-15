#include <gtest/gtest.h>
#include <ips4o.hpp>
#include <lz77_sss/data_structures/sample_index/sample_index.hpp>

using lce_r_t = alx::lce::lce_naive_wordwise_xor<char>;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<uint32_t> avg_sample_rate_distrib(1, 10);
std::uniform_int_distribution<uint32_t> pattern_length_distrib(1, 1000);
std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);

std::string input;
uint32_t avg_sample_rate;
std::vector<uint32_t> sampling;
sample_index<>::query_context query;
uint32_t pattern_pos;
uint32_t pattern_length;
std::vector<uint32_t> occurrences;
std::vector<uint32_t> correct_occurrences;

template <direction dir>
void test_query(sample_index<>& index) {
    std::uniform_int_distribution<uint32_t> pattern_pos_distrib(0, input.size() - 1);
    pattern_pos = pattern_pos_distrib(gen);
    query = index.query();
    pattern_length = std::min<uint32_t>(
        dir == LEFT ? pattern_pos + 1 : (input.size() - pattern_pos),
        pattern_length_distrib(gen));
    std::uniform_int_distribution<uint32_t> step_size_distrib(1, 3);
    uint32_t current_pattern_length = 0;
    bool occurs = true;

    while (occurs && current_pattern_length < pattern_length) {
        current_pattern_length = std::min<uint32_t>(
            pattern_length, current_pattern_length + step_size_distrib(gen));
        occurs = index.extend<dir>(query, pattern_pos, current_pattern_length);
    }

    if (occurs) {
        index.locate<dir>(query, occurrences);
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
    correct_occurrences.clear();
    occurrences.clear();
}

TEST(test_sample_index, fuzzy_test) {
    auto start_time = now();

    while (time_diff_min(start_time,now()) < 60) {
        // generate a random string
        input = random_repetitive_string(1, 100000);

        // choose a random average sample rate
        avg_sample_rate = avg_sample_rate_distrib(gen);
        std::uniform_int_distribution<uint32_t> sample_distance_distrib(1, 2 * avg_sample_rate);

        // compute a random sampling of text positions
        sampling.emplace_back(std::min<uint32_t>(input.size() - 1, sample_distance_distrib(gen)));
        while (sampling.back() + 2 * avg_sample_rate < input.size()) {
            sampling.emplace_back(sampling.back() + sample_distance_distrib(gen));
        }

        // build the sample-index
        sample_index<> index;
        index.build(input, sampling, lce_r_t(input), true, omp_get_max_threads());

        // perform random queries and check their correctness
        for (uint32_t i = 0; i < 1000; i++) {
            if (prob_distrib(gen) < 0.5) {
                test_query<LEFT>(index);
            } else {
                test_query<RIGHT>(index);
            }
        }
        
        sampling.clear();
        input.clear();
    }
}