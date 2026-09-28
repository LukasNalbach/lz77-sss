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

#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <lz77_sss/misc/utils.hpp>

enum class repetitiveness_kind : uint8_t {
    periodic,
    versioned,
    indel,
    block_move,
    markov,
    lz,
    fibonacci,
    count
};

inline std::string_view repetitiveness_kind_name(repetitiveness_kind kind)
{
    switch (kind) {
        case repetitiveness_kind::periodic: return "periodic";
        case repetitiveness_kind::versioned: return "versioned";
        case repetitiveness_kind::indel: return "indel";
        case repetitiveness_kind::block_move: return "blockmove";
        case repetitiveness_kind::markov: return "markov";
        case repetitiveness_kind::lz: return "lz";
        default: return "fibonacci";
    }
}

struct repetitive_params {
    repetitiveness_kind kind = repetitiveness_kind::versioned;
    uint64_t size = 1 << 20;
    uint64_t base_length = 4096;
    uint64_t block_length = 256;
    double mutation_rate = 0.001;
    uint64_t seed = 0;
};

template <typename inp_t>
static inp_t generate_repetitive_input(
    repetitive_params params,
    typename inp_t::value_type min_sym = std::numeric_limits<typename inp_t::value_type>::min(),
    typename inp_t::value_type max_sym = std::numeric_limits<typename inp_t::value_type>::max())
{
    using sym_t = typename inp_t::value_type;
    using draw_t = std::conditional_t<sizeof(sym_t) == 1,
        std::conditional_t<std::is_signed_v<sym_t>, int32_t, uint32_t>, sym_t>;

    std::random_device device;
    std::mt19937_64 random(params.seed == 0 ? device() : params.seed);
    std::uniform_int_distribution<draw_t> symbols(min_sym, max_sym);
    std::uniform_real_distribution<double> chance(0.0, 1.0);

    params.size = std::max<uint64_t>(1, params.size);
    params.base_length = std::clamp<uint64_t>(params.base_length, 1, params.size);
    params.block_length = std::clamp<uint64_t>(params.block_length, 1, params.base_length);

    inp_t input;
    input.reserve(params.size);
    auto draw = [&]() { return sym_t(symbols(random)); };

    if (params.kind == repetitiveness_kind::fibonacci) {
        inp_t previous, current;
        previous.push_back(sym_t(min_sym));
        current.push_back(sym_t(max_sym == min_sym ? min_sym : min_sym + 1));

        while (current.size() < params.size) {
            inp_t next = current;
            next.insert(next.end(), previous.begin(), previous.end());
            previous = std::move(current);
            current = std::move(next);
        }

        current.resize(params.size);
        return current;
    }

    if (params.kind == repetitiveness_kind::markov) {
        uint64_t least = std::max<uint64_t>(1, params.base_length / 4);

        while (input.size() < params.size) {
            if (input.size() > params.base_length && chance(random) < 0.5) {
                std::uniform_int_distribution<uint64_t> from(0, input.size() - params.base_length - 1);
                std::uniform_int_distribution<uint64_t> length(least, params.base_length);
                uint64_t start = from(random);
                uint64_t take = std::min(length(random), params.size - input.size());
                input.insert(input.end(), input.begin() + start, input.begin() + start + take);
            } else {
                for (uint64_t i = 0; i < least && input.size() < params.size; i++)
                    input.push_back(draw());
            }
        }

        input.resize(params.size);
        return input;
    }

    if (params.kind == repetitiveness_kind::lz) {
        input.push_back(draw());
        std::uniform_int_distribution<uint64_t> length(1, params.base_length);

        while (input.size() < params.size) {
            if (chance(random) < 0.98) {
                std::uniform_int_distribution<uint64_t> from(0, input.size() - 1);
                uint64_t start = from(random);
                uint64_t take = std::min({ length(random), input.size() - start,
                    params.size - input.size() });
                input.insert(input.end(), input.begin() + start, input.begin() + start + take);
            } else {
                input.push_back(draw());
            }
        }

        input.resize(params.size);
        return input;
    }

    inp_t base;
    base.reserve(params.base_length);
    for (uint64_t i = 0; i < params.base_length; i++) base.push_back(draw());

    while (input.size() < params.size) {
        switch (params.kind) {
            case repetitiveness_kind::periodic:
                input.insert(input.end(), base.begin(),
                    base.begin() + std::min<uint64_t>(base.size(), params.size - input.size()));
                break;

            case repetitiveness_kind::versioned:
                for (uint64_t i = 0; i < base.size() && input.size() < params.size; i++)
                    input.push_back(chance(random) < params.mutation_rate ? draw() : base[i]);
                break;

            case repetitiveness_kind::indel:
                for (uint64_t i = 0; i < base.size() && input.size() < params.size; i++) {
                    double roll = chance(random);
                    if (roll < params.mutation_rate / 2) continue;
                    if (roll < params.mutation_rate) input.push_back(draw());
                    input.push_back(base[i]);
                }
                break;

            default: {
                uint64_t blocks = std::max<uint64_t>(1, base.size() / params.block_length);
                std::vector<uint64_t> order(blocks);
                for (uint64_t i = 0; i < blocks; i++) order[i] = i;
                std::shuffle(order.begin(), order.end(), random);

                for (uint64_t block : order) {
                    if (input.size() >= params.size) break;
                    uint64_t start = block * params.block_length;
                    uint64_t take = std::min({ params.block_length, base.size() - start,
                        params.size - input.size() });
                    input.insert(input.end(), base.begin() + start, base.begin() + start + take);
                }
                break;
            }
        }
    }

    input.resize(params.size);
    return input;
}

template <typename rng_t>
static repetitive_params random_repetitive_params(uint64_t size, rng_t& random)
{
    std::uniform_int_distribution<uint32_t> kinds(0, uint32_t(repetitiveness_kind::count) - 1);
    std::uniform_real_distribution<double> chance(0.0, 1.0);
    repetitive_params params;

    params.kind = repetitiveness_kind(kinds(random));
    params.size = size;

    std::uniform_int_distribution<uint64_t> bases(1, std::max<uint64_t>(1, size));
    params.base_length = std::max<uint64_t>(1, bases(random) / (1 + kinds(random)));
    params.block_length = std::max<uint64_t>(1, params.base_length / (1 + kinds(random)));
    params.mutation_rate = chance(random) * 0.05;
    params.seed = std::uniform_int_distribution<uint64_t>(1, ~uint64_t(0))(random);

    return params;
}

template <typename inp_t>
static inp_t random_repetitive_input(
    uint64_t min_size, uint64_t max_size,
    typename inp_t::value_type min_sym = std::numeric_limits<typename inp_t::value_type>::min(),
    typename inp_t::value_type max_sym = std::numeric_limits<typename inp_t::value_type>::max(),
    uint64_t seed = 0)
{
    std::random_device device;
    std::mt19937_64 random(seed == 0 ? device() : seed);
    uint64_t size = random_log_uniform_size(min_size, max_size, random);

    inp_t input = generate_repetitive_input<inp_t>(
        random_repetitive_params(size, random), min_sym, max_sym);

    if constexpr (std::is_same_v<inp_t, std::string>) {
        no_init_resize_with_excess(input, input.size(), 4 * 4096);
    }

    return input;
}
