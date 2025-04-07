#include <gtest/gtest.h>
#include <ips4o.hpp>
#include <lz77_sss/data_structures/dynamic_range/dynamic_square_grid.hpp>
#include <lz77_sss/data_structures/dynamic_range/semi_dynamic_square_grid.hpp>

using point_t = dynamic_range<>::point_t;

struct query {
    uint32_t x1, x2;
    uint32_t y1, y2;
    bool result;
};

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<uint32_t> input_size_distrib(1, 10000);
std::uniform_int_distribution<uint32_t> range_size_distrib(1, 100000);

uint32_t input_size;
std::vector<point_t> input;
uint32_t pos_max;
std::vector<query> queries;
uint32_t num_queries;

template <typename range_ds_t>
void test()
{
    // choose a random input length
    input_size = input_size_distrib(gen);

    // choose a random range to pick points from
    pos_max = range_size_distrib(gen);
    std::uniform_int_distribution<uint32_t> pos_distrib(0, pos_max - 1);

    // generate a random set of input points
    for (uint32_t i = 0; i < input_size; i++) {
        input.emplace_back(point_t {
            .x = pos_distrib(gen),
            .y = pos_distrib(gen) });
    }

    // remove duplicate points
    input.erase(std::unique(input.begin(), input.end(),
        [](point_t& p1, point_t& p2) {
            return p1.x == p2.x && p1.y == p2.y;
        }),
        input.end());

    // generate random queries
    std::vector<query> queries;
    num_queries = input.size();
    for (uint32_t i = 0; i < num_queries; i++) {
        query q {
            .x1 = pos_distrib(gen),
            .x2 = pos_distrib(gen),
            .y1 = pos_distrib(gen),
            .y2 = pos_distrib(gen),
            .result = false
        };

        if (q.x1 > q.x2) std::swap(q.x1, q.x2);
        if (q.y1 > q.y2) std::swap(q.y1, q.y2);

        for (uint32_t j = 0; j < i; j++) {
            const point_t& p = input[j];
            if (q.x1 <= p.x && p.x <= q.x2 &&
                q.y1 <= p.y && p.y <= q.y2) {
                q.result = true;
                break;
            }
        }

        queries.emplace_back(q);
    }

    // build the range data structure
    range_ds_t ds(input, pos_max);

    // verify that all queries are answered correctly
    for (uint32_t i = 0; i < num_queries; i++) {
        const query& q = queries[i];
        auto [p, result] = ds.point_in_range(
            q.x1, q.x2, q.y1, q.y2);
        EXPECT_EQ(result, q.result);
        EXPECT_TRUE(!result ||
            (q.x1 <= p.x && p.x <= q.x2 &&
            q.y1 <= p.y && p.y <= q.y2));
        ds.insert(input[i]);
    }

    input.clear();
    queries.clear();
}

TEST(test_dynamic_range, fuzzy_test)
{
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        switch (std::rand() % 2) {
            case 0: test<dynamic_square_grid<>>(); break;
            case 1: test<semi_dynamic_square_grid<>>(); break;
        }
    }
}