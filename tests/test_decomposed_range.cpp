#include <gtest/gtest.h>
#include <ips4o.hpp>
#include <lz77_sss/data_structures/dynamic_range/dynamic_square_grid.hpp>
#include <lz77_sss/data_structures/dynamic_range/semi_dynamic_square_grid.hpp>
#include <lz77_sss/data_structures/sample_index/sample_index.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_kd_tree.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_square_grid.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_striped_square.hpp>

using point_t = static_weighted_range<>::point_t;
using lce_r_t = lce::ds::lce_naive_wordwise_xor<char>;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<uint32_t> input_size_distrib(10000, 100000);
std::uniform_int_distribution<uint32_t> avg_sample_rate_distrib(1, 10);

struct query {
    char chr;
    uint32_t x1, x2;
    uint32_t y1, y2;
    uint32_t weight;
    bool result;
};

template <template <typename> typename range_ds_t>
void test()
{
    using point_t = typename range_ds_t<uint32_t>::point_t;

    // choose a random input length
    uint32_t input_size = input_size_distrib(gen);

    // generate a random string
    char* input = nullptr;
    random_repetitive_string(input, input_size, input_size);

    // choose a random average sample rate
    uint32_t avg_sample_rate = avg_sample_rate_distrib(gen);
    std::uniform_int_distribution<uint32_t> sample_distance_distrib(1, 2 * avg_sample_rate);

    // compute a random sampling of text positions
    std::vector<uint32_t> sampling;
    sampling.emplace_back(std::min<uint32_t>(input_size - 1, sample_distance_distrib(gen)));
    while (sampling.back() + 2 * avg_sample_rate < input_size) {
        sampling.emplace_back(sampling.back() + sample_distance_distrib(gen));
    }
    uint32_t num_samples = sampling.size();

    // build a sample index (SA_S and PA_S)
    sample_index<> index;
    index.build(input, input_size, sampling, lce_r_t(input, input_size), false);

    // build the points-array
    std::vector<point_t> points;
    points.reserve(num_samples);

    for (uint32_t i = 0; i < num_samples; i++) {
        if constexpr (range_ds_t<uint32_t>::is_static()) {
            points.emplace_back(point_t { .weight = i });
        } else {
            points.emplace_back(point_t { });
        }
    }

    for (uint32_t i = 0; i < num_samples; i++) {
        points[index.pa_s(i)].x = i;
        points[index.sa_s(i)].y = i;
    }

    // generate random queries
    std::vector<query> queries;
    std::uniform_int_distribution<uint8_t> uchar_distrib(0, 255);
    std::array<std::uniform_int_distribution<uint32_t>, 256> query_range_distrib;
    std::array<uint32_t, 257> c_array = { 0 };
    std::vector<uint8_t> used_uchars;

    for (uint32_t sample : sampling) {
        c_array[char_to_uchar(input[sample])]++;
    }

    for (uint16_t c = 1; c < 256; c++) c_array[c] += c_array[c - 1];
    for (uint16_t c = 256; c > 0; c--) c_array[c] = c_array[c - 1];
    c_array[0] = 0;

    for (uint16_t c = 0; c < 256; c++) {
        if (c_array[c] != c_array[c + 1]) {
            used_uchars.emplace_back(c);
        }
    }

    std::uniform_int_distribution<uint8_t>
        uchar_idx_distrib(0, used_uchars.size() - 1);

    for (uint16_t c = 0; c < 256; c++) {
        query_range_distrib[c] = std::uniform_int_distribution<uint32_t>(
            c_array[c], c_array[c + 1] - 1);
    }

    for (uint32_t i = 0; i < num_samples; i++) {
        uint8_t uchar = used_uchars[uchar_idx_distrib(gen)];

        query q {
            .chr = uchar_to_char<char>(uchar),
            .x1 = query_range_distrib[uchar](gen),
            .x2 = query_range_distrib[uchar](gen),
            .y1 = query_range_distrib[uchar](gen),
            .y2 = query_range_distrib[uchar](gen),
            .weight = i,
            .result = false
        };

        if (q.x1 > q.x2) std::swap(q.x1, q.x2);
        if (q.y1 > q.y2) std::swap(q.y1, q.y2);

        for (uint32_t j = 0; j < i; j++) {
            const point_t& p = points[j];
            if (q.x1 <= p.x && p.x <= q.x2 &&
                q.y1 <= p.y && p.y <= q.y2
            ) {
                q.result = true;
                break;
            }
        }

        queries.emplace_back(q);
    }

    // build the range data structure
    range_ds_t<uint32_t> ds(input, sampling, points);
    if constexpr (range_ds_t<uint32_t>::is_static()) points.clear();

    // verify that all queries are answered correctly
    for (uint32_t i = 0; i < num_samples; i++) {
        const query& q = queries[i];
        point_t p;
        bool result;

        if constexpr (range_ds_t<uint32_t>::is_static()) {
            std::tie(p, result) = ds.lighter_point_in_range(
                q.chr, q.weight, q.x1, q.x2, q.y1, q.y2);
        } else {
            std::tie(p, result) = ds.point_in_range(
                q.chr, q.x1, q.x2, q.y1, q.y2);
            ds.insert(input[sampling[i]], points[i]);
        }

        EXPECT_EQ(result, q.result);

        if (result) {
            EXPECT_TRUE(
                q.x1 <= p.x && p.x <= q.x2 &&
                q.y1 <= p.y && p.y <= q.y2);

            if constexpr (range_ds_t<uint32_t>::is_static()) {
                EXPECT_TRUE(p.weight < q.weight);
            }
        }
    }

    queries.clear();
    sampling.clear();
}

TEST(test_static_weighted_range, fuzzy_test)
{
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        std::string input;
        std::vector<uint32_t> samples;

        switch (std::rand() % 5) {
            case 0: test<decomposed_static_weighted_kd_tree>(); break;
            case 1: test<decomposed_static_weighted_square_grid>(); break;
            case 2: test<decomposed_static_weighted_striped_square>(); break;
            case 3: test<decomposed_dynamic_square_grid>(); break;
            case 4: test<decomposed_semi_dynamic_square_grid>(); break;
        }
    }
}