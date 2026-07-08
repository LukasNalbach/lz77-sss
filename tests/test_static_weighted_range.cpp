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
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_kd_tree.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_square_grid.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_striped_square.hpp>

#include "test-progress.hpp"

using point_t = static_weighted_range<>::point_t;

struct query {
    uint32_t x1, x2;
    uint32_t y1, y2;
    uint32_t weight;
    bool result;
};

thread_local std::mt19937 gen(std::random_device{}());

template <typename range_ds_t>
void test()
{
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);

    // choose a random input length
    uint32_t input_size = random_log_uniform_size(1, 10000, gen);
    std::vector<point_t> input;
    std::uniform_int_distribution<uint32_t> distrib(0, input_size - 1);

    // generate a random set of input points with random weights,
    // s.t. the sets of x-coordinates, y-coordinates and weights
    // are permutations of [0, input_size - 1]
    input.resize(input_size);
    for (uint32_t i = 0; i < input_size; i++) input[i].x = i;
    std::shuffle(input.begin(), input.end(), gen);
    for (uint32_t i = 0; i < input_size; i++) input[i].y = i;
    std::shuffle(input.begin(), input.end(), gen);
    for (uint32_t i = 0; i < input_size; i++) input[i].weight = i;
    std::shuffle(input.begin(), input.end(), gen);

    // generate random queries
    std::vector<query> queries(1000);
    for (query& q : queries) {
        q.x1 = distrib(gen);
        q.x2 = distrib(gen);
        q.y1 = distrib(gen);
        q.y2 = distrib(gen);
        q.weight = distrib(gen);
        q.result = false;
        if (q.x1 > q.x2) std::swap(q.x1, q.x2);
        if (q.y1 > q.y2) std::swap(q.y1, q.y2);
    }

    #pragma omp parallel for num_threads(num_threads)
    for (int64_t k = 0; k < (int64_t) queries.size(); k++) {
        query& q = queries[k];
        for (const point_t& p : input) {
            if (p.weight < q.weight &&
                q.x1 <= p.x && p.x <= q.x2 &&
                q.y1 <= p.y && p.y <= q.y2
            ) {
                q.result = true;
                break;
            }
        }
    }

    // build the range data structure
    range_ds_t ds(input, input_size, num_threads);

    // verify that all queries are answered correctly
    #pragma omp parallel for num_threads(num_threads)
    for (int64_t k = 0; k < (int64_t) queries.size(); k++) {
        query& q = queries[k];
        auto [p, result] = ds.lighter_point_in_range(
            q.weight, q.x1, q.x2, q.y1, q.y2);
        EXPECT_EQ(result, q.result);
        EXPECT_TRUE(!result ||
            (p.weight < q.weight &&
            q.x1 <= p.x && p.x <= q.x2 &&
            q.y1 <= p.y && p.y <= q.y2));
    }
}

TEST(test_static_weighted_range, static_weighted_kd_tree)
{
    run_fuzz("static-weighted-range", {
        { "static-weighted-kd-tree", [](uint64_t) { test<static_weighted_kd_tree<>>(); }, false },
    }, fuzz_iterations(6000));
}

TEST(test_static_weighted_range, static_weighted_square_grid)
{
    run_fuzz("static-weighted-range", {
        { "static-weighted-square-grid", [](uint64_t) { test<static_weighted_square_grid<>>(); }, false },
    }, fuzz_iterations(6000));
}

TEST(test_static_weighted_range, static_weighted_striped_square)
{
    run_fuzz("static-weighted-range", {
        { "static-weighted-striped-square", [](uint64_t) { test<static_weighted_striped_square<>>(); }, false },
    }, fuzz_iterations(6000));
}