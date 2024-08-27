#include <gtest/gtest.h>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_kd_tree.hpp>

using point_t = static_weighted_kd_tree<>::point_t;

struct query {
    uint32_t x1, x2;
    uint32_t y1, y2;
    uint32_t weight;
    bool result;
};

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<uint32_t> input_size_distrib(1, 10000);
std::uniform_int_distribution<uint32_t> range_size_distrib(1, 100000);
std::uniform_int_distribution<uint32_t> weight_distrib(0, 10000);

uint32_t input_size;
std::vector<point_t> input;
uint32_t x_min, x_max, y_min, y_max;
std::vector<query> queries;

TEST(test_static_weighted_kd_tree, fuzzy_test) {
    auto start_time = now();

    while (time_diff_min(start_time,now()) < 60) {
        // choose a random input length
        input_size = input_size_distrib(gen);

        // choose a random range to pick points from
        x_min = range_size_distrib(gen);
        x_max = range_size_distrib(gen);
        y_min = range_size_distrib(gen);
        x_max = range_size_distrib(gen);
        if (x_min > x_max) std::swap(x_min, x_max);
        if (y_min > y_max) std::swap(y_min, y_max);
        std::uniform_int_distribution<uint32_t> x_distrib(x_min, x_max);
        std::uniform_int_distribution<uint32_t> y_distrib(y_min, y_max);

        // generate a random set of input points with random weights
        for (uint32_t i = 0; i < input_size; i++) {
            input.emplace_back(point_t{
                .x = x_distrib(gen),
                .y = y_distrib(gen),
                .weight = weight_distrib(gen)
            });
        }

        // remove duplicate points
        input.erase(std::unique(input.begin(), input.end(),
            [](point_t& p1, point_t& p2){
                return p1.x == p2.x && p1.y == p2.y;
        }), input.end());

        // generate random queries
        std::vector<query> queries;
        for (uint32_t i = 0; i < 1000; i++) {
            query q {
                .x1 = x_distrib(gen),
                .x2 = x_distrib(gen),
                .y1 = y_distrib(gen),
                .y2 = y_distrib(gen),
                .weight = weight_distrib(gen),
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

        // build the tree
        static_weighted_kd_tree<> tree(std::move(input));

        // verify that all queries are answered correctly
        for (query& q : queries) {
            auto [p, result] = tree.lighter_point_in_range(
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
}