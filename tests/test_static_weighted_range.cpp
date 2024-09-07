#include <gtest/gtest.h>
#include <ips4o.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_kd_tree.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_square_grid.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_striped_square.hpp>

using point_t = static_weighted_range<>::point_t;

struct query {
    uint32_t x1, x2;
    uint32_t y1, y2;
    uint32_t weight;
    bool result;
};

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<uint32_t> input_size_distrib(1, 10000);

uint32_t input_size;
std::vector<point_t> input;
std::vector<query> queries;

template <typename range_ds_t>
void test() {
    // choose a random input length
    input_size = input_size_distrib(gen);
    std::uniform_int_distribution<uint32_t> distrib(0, input_size - 1);

    // generate a random set of input points with random weights,
    // s.t. the sets of x-coordinates, y-coordinates and weights
    // are permutations of [0, input_size - 1]
    input.resize(input_size);
    for (uint32_t i = 0; i < input_size; i++) input[i].x = i;
    std::random_shuffle(input.begin(), input.end());
    for (uint32_t i = 0; i < input_size; i++) input[i].y = i;
    std::random_shuffle(input.begin(), input.end());
    for (uint32_t i = 0; i < input_size; i++) input[i].weight = i;
    std::random_shuffle(input.begin(), input.end());

    // generate random queries
    std::vector<query> queries;
    for (uint32_t i = 0; i < 1000; i++) {
        query q {
            .x1 = distrib(gen),
            .x2 = distrib(gen),
            .y1 = distrib(gen),
            .y2 = distrib(gen),
            .weight = distrib(gen),
            .result = false
        };

        if (q.x1 > q.x2) std::swap(q.x1, q.x2);
        if (q.y1 > q.y2) std::swap(q.y1, q.y2);

        for (point_t& p : input) {
            if (p.weight < q.weight &&
                q.x1 <= p.x && p.x <= q.x2 &&
                q.y1 <= p.y && p.y <= q.y2
            ) {
                q.result = true;
                break;
            }
        }

        queries.emplace_back(q);
    }

    // build the range data structure
    range_ds_t ds(input, input_size);

    // verify that all queries are answered correctly
    for (query& q : queries) {
        auto [p, result] = ds.lighter_point_in_range(
            q.weight, q.x1, q.x2, q.y1, q.y2);
        EXPECT_EQ(result, q.result);
        EXPECT_TRUE(!result || (
            p.weight < q.weight &&
            q.x1 <= p.x && p.x <= q.x2 &&
            q.y1 <= p.y && p.y <= q.y2
        ));
    }

    queries.clear();
}

TEST(test_static_weighted_range, fuzzy_test) {
    auto start_time = now();

    while (time_diff_min(start_time,now()) < 60) {
        switch (std::rand() % 3) {
            case 0: test<static_weighted_kd_tree<>>(); break;
            case 1: test<static_weighted_square_grid<>>(); break;
            case 2: test<static_weighted_striped_square<>>(); break;
        }
    }
}