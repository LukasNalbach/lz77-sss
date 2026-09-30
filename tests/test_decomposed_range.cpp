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
#include <map>
#include <ips4o.hpp>
#include <lz77_sss/data_structures/range/range.hpp>
#include <lz77_sss/data_structures/sample_index/sample_index.hpp>
#include <lz77_sss/misc/repetitive_input.hpp>

#include "test-progress.hpp"

using point_t = range_ds::point_t;

thread_local std::mt19937 gen(std::random_device{}());
thread_local std::uniform_int_distribution<uint32_t> avg_sample_rate_distrib(1, 10);

struct query {
    uint64_t chr;
    uint32_t x1, x2;
    uint32_t y1, y2;
    uint32_t weight;
    bool result;
};

template <typename input_t>
void test(const range_ds_kind& kind)
{
    using sym_t = typename input_t::value_type;
    static constexpr bool byte_input = sizeof(sym_t) == 1;
    using text_t = lce::text::direct_text<sym_t>;
    using lce_r_t = lce::ds::lce_naive_wordwise_xor<sym_t>;
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);

    // choose a random input length
    uint32_t input_size = random_log_uniform_size(1, 10000, gen);

    // generate a random string
    input_t input;

    if constexpr (byte_input) {
        input = random_repetitive_input<input_t>(input_size, input_size);
    } else {
        const uint32_t sigma = random_log_uniform_size(1, uint64_t { 1 } << 16, gen);
        input = random_repetitive_input<input_t>(input_size, input_size, sym_t(0), sym_t(sigma - 1));
    }

    auto symbol_at = [&](uint64_t i) -> uint64_t {
        if constexpr (byte_input) return char_to_uchar(input[i]);
        else return input[i];
    };

    // choose a random average sample rate
    uint32_t avg_sample_rate = avg_sample_rate_distrib(gen);
    std::uniform_int_distribution<uint32_t> sample_distance_distrib(1, 2 * avg_sample_rate);

    // compute a random sampling of text positions
    std::vector<uint40_t> sampling;
    sampling.emplace_back(std::min<uint64_t>(input_size - 1, sample_distance_distrib(gen)));
    while (uint64_t(sampling.back()) + 2 * avg_sample_rate < input_size) {
        sampling.emplace_back(uint64_t(sampling.back()) + sample_distance_distrib(gen));
    }
    uint32_t num_samples = sampling.size();

    // build a sample index (SA_S and PA_S)
    text_t text;

    if constexpr (byte_input) {
        text = text_t(input.data(), input_size);
    } else {
        text = text_t(input.data(), input_size, *std::max_element(input.begin(), input.end()) + uint64_t(1));
    }

    sample_index<text_t, lce_r_t> index;
    const lce_r_t lce_r(input);
    index.build(text, input_size, sampling, lce_r, interval_samples::skip, 32, num_threads);

    // build the points-array
    std::vector<point_t> points;
    points.reserve(num_samples);

    for (uint32_t i = 0; i < num_samples; i++) {
        if (kind.is_static()) {
            points.emplace_back(point_t { .x = 0, .y = 0, .weight = i });
        } else {
            points.emplace_back(point_t { });
        }
    }

    for (uint32_t i = 0; i < num_samples; i++) {
        points[index.pa(i)].x = i;
        points[index.sa(i)].y = i;
    }

    // generate random queries
    std::vector<query> queries;
    std::map<uint64_t, uint32_t> block_beg;
    std::map<uint64_t, uint32_t> block_len;
    for (uint32_t sample : sampling) block_len[symbol_at(sample)]++;
    std::vector<uint64_t> used_symbols;
    uint32_t sum = 0;

    for (auto [c, len] : block_len) {
        block_beg[c] = sum;
        sum += len;
        used_symbols.emplace_back(c);
    }

    std::uniform_int_distribution<uint64_t> symbol_idx_distrib(0, used_symbols.size() - 1);

    for (uint32_t i = 0; i < num_samples; i++) {
        const uint64_t c = used_symbols[symbol_idx_distrib(gen)];
        std::uniform_int_distribution<uint32_t> query_range_distrib(block_beg[c], block_beg[c] + block_len[c] - 1);

        query q {
            .chr = c,
            .x1 = query_range_distrib(gen),
            .x2 = query_range_distrib(gen),
            .y1 = query_range_distrib(gen),
            .y2 = query_range_distrib(gen),
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
    range_ds* ds = make_range_ds(kind, input.data(), sampling, points, num_threads);

    // verify that all queries are answered correctly
    for (uint32_t i = 0; i < num_samples; i++) {
        const query& q = queries[i];
        point_t p;
        bool result;

        if (kind.is_static()) {
            std::tie(p, result) = ds->lighter_point_in_range(
                q.chr, q.weight, q.x1, q.x2, q.y1, q.y2);
        } else {
            std::tie(p, result) = ds->point_in_range(
                q.chr, q.x1, q.x2, q.y1, q.y2);
            ds->insert(symbol_at(sampling[i]), points[i]);
        }

        EXPECT_EQ(result, q.result);

        if (result) {
            EXPECT_TRUE(
                q.x1 <= p.x && p.x <= q.x2 &&
                q.y1 <= p.y && p.y <= q.y2);

            if (kind.is_static()) {
                EXPECT_TRUE(uint64_t(p.weight) < q.weight);
            }
        }
    }

    delete ds;
}

TEST(test_decomposed_range, decomposed_static_weighted_kd_tree)
{
    run_fuzz("decomposed-range", {
        { "decomposed-static-weighted-kd-tree", [](uint64_t) { test<std::string>(range_ds_kind { .type = range_ds_type::swkdt, .decomposed = true }); }, false },
    }, fuzz_iterations(3000));
}

TEST(test_decomposed_range, decomposed_static_weighted_square_grid)
{
    run_fuzz("decomposed-range", {
        { "decomposed-static-weighted-square-grid", [](uint64_t) { test<std::string>(range_ds_kind { .type = range_ds_type::swsg, .decomposed = true }); }, false },
    }, fuzz_iterations(3000));
}

TEST(test_decomposed_range, decomposed_semi_dynamic_square_grid)
{
    run_fuzz("decomposed-range", {
        { "decomposed-semi-dynamic-square-grid", [](uint64_t) { test<std::string>(range_ds_kind { .type = range_ds_type::sdsg, .decomposed = true }); }, false },
    }, fuzz_iterations(3000));
}

TEST(test_decomposed_range, grouped_static_weighted_kd_tree)
{
    run_fuzz("decomposed-range", {
        { "grouped-static-weighted-kd-tree", [](uint64_t) { test<std::vector<uint32_t>>(range_ds_kind { .type = range_ds_type::swkdt, .decomposed = true }); }, false },
    }, fuzz_iterations(3000));
}

TEST(test_decomposed_range, grouped_static_weighted_square_grid)
{
    run_fuzz("decomposed-range", {
        { "grouped-static-weighted-square-grid", [](uint64_t) { test<std::vector<uint32_t>>(range_ds_kind { .type = range_ds_type::swsg, .decomposed = true }); }, false },
    }, fuzz_iterations(3000));
}

TEST(test_decomposed_range, grouped_semi_dynamic_square_grid)
{
    run_fuzz("decomposed-range", {
        { "grouped-semi-dynamic-square-grid", [](uint64_t) { test<std::vector<uint32_t>>(range_ds_kind { .type = range_ds_type::sdsg, .decomposed = true }); }, false },
    }, fuzz_iterations(3000));
}
