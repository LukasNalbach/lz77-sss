/**
 * part of LukasNalbach/lz77-sss
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <gtest/gtest.h>
#include <ips4o.hpp>
#include <lz77_sss/data_structures/dynamic_range/dynamic_square_grid.hpp>
#include <lz77_sss/data_structures/dynamic_range/semi_dynamic_square_grid.hpp>

#include "test-progress.hpp"

using point_t = dynamic_range<>::point_t;

struct query {
    uint32_t x1, x2;
    uint32_t y1, y2;
    bool result;
};

thread_local std::mt19937 gen(std::random_device{}());
thread_local std::uniform_int_distribution<uint32_t> range_size_distrib(1, 100000);

template <typename range_ds_t>
void test()
{
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);

    // choose a random input length
    uint32_t input_size = random_log_uniform_size(1, 10000, gen);
    std::vector<point_t> input;

    // choose a random range to pick points from
    uint32_t pos_max = range_size_distrib(gen);
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
    uint32_t num_queries = input.size();
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
    range_ds_t ds(input, pos_max, num_threads);

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
}

TEST(test_dynamic_range, dynamic_square_grid)
{
    run_fuzz("dynamic-range", {
        { "dynamic-square-grid", [](uint64_t) { test<dynamic_square_grid<>>(); }, false },
    }, fuzz_iterations(6000));
}

TEST(test_dynamic_range, semi_dynamic_square_grid)
{
    run_fuzz("dynamic-range", {
        { "semi-dynamic-square-grid", [](uint64_t) { test<semi_dynamic_square_grid<>>(); }, false },
    }, fuzz_iterations(6000));
}