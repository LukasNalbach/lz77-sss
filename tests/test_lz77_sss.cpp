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
#include <unordered_set>
#include <lz77_sss/inst.hpp>
#include <lz77_sss/misc/fasta_interleaver.hpp>
#include <lz77_sss/misc/file_decoder.hpp>
#include <lz77_sss/misc/text_reader.hpp>
#include <lz77_sss/misc/repetitive_input.hpp>
#include <lz77/parallel_lpf.hpp>

#include "test-progress.hpp"

std::random_device rd;
std::mt19937 gen(rd());

using factor_out = std::function<void(lz77_sss::factor)>;

void run_approximate(std::string& in, const factor_out& out, uint16_t num_threads) {
    lz77_sss::factorize_approximate<>(in.data(), in.size(), out, { .num_threads = num_threads, .fact_mode = lz77_sss::auto_gaps });
}

void run_approximate_exact_gaps(std::string& in, const factor_out& out, uint16_t num_threads) {
    lz77_sss::factorize_approximate<>(in.data(), in.size(), out, { .num_threads = num_threads, .fact_mode = lz77_sss::exact_gaps });
}

void run_exact_with_interval_samples(std::string& in, const factor_out& out, uint16_t num_threads) {
    switch (std::rand() % 2) {
        case 0: lz77_sss::factorize_exact<>(in.data(), in.size(), out, { .num_threads = num_threads, .fact_mode = lz77_sss::auto_gaps, .transf_mode = lz77_sss::with_interval_samples, .range_ds = { .type = range_ds_type::sdsg, .decomposed = true }, .exact_alg = lz77_sss::sss_based }); break;
        case 1: lz77_sss::factorize_exact<>(in.data(), in.size(), out, { .num_threads = num_threads, .fact_mode = lz77_sss::auto_gaps, .transf_mode = lz77_sss::with_interval_samples, .range_ds = { .type = range_ds_type::swkdt, .decomposed = true }, .exact_alg = lz77_sss::sss_based }); break;
    }
}

void run_exact_without_interval_samples(std::string& in, const factor_out& out, uint16_t num_threads) {
    switch (std::rand() % 2) {
        case 0: lz77_sss::factorize_exact<>(in.data(), in.size(), out, { .num_threads = num_threads, .fact_mode = lz77_sss::auto_gaps, .transf_mode = lz77_sss::without_interval_samples, .range_ds = { .type = range_ds_type::sdsg, .decomposed = true }, .exact_alg = lz77_sss::sss_based }); break;
        case 1: lz77_sss::factorize_exact<>(in.data(), in.size(), out, { .num_threads = num_threads, .fact_mode = lz77_sss::auto_gaps, .transf_mode = lz77_sss::without_interval_samples, .range_ds = { .type = range_ds_type::swkdt, .decomposed = true }, .exact_alg = lz77_sss::sss_based }); break;
    }
}

void run_exact_sa(std::string& in, const factor_out& out, uint16_t num_threads) {
    lz77_sss::factorize_exact<>(in.data(), in.size(), out, { .num_threads = num_threads, .exact_alg = lz77_sss::sa_based });
}

void run_exact_auto(std::string& in, const factor_out& out, uint16_t num_threads) {
    lz77_sss::factorize_exact<>(in.data(), in.size(), out, { .num_threads = num_threads, .exact_alg = lz77_sss::auto_select });
}

std::vector<uint64_t> reference_phrase_lengths(const std::string& input) {
    std::vector<uint64_t> lengths;
    lz77::emit_function emit = [&](lz77::factor f) { lengths.push_back(f.text_len()); };
    lz77::parallel_lpf_factorizer().factorize(input.begin(), input.end(), emit, emit, 1);
    return lengths;
}

void test_roundtrip(std::string input, const std::function<void(std::string&, const factor_out&, uint16_t)>& run, bool exact = false) {
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    std::vector<lz77_sss::factor> factorization;
    run(input, [&](lz77_sss::factor f) { factorization.emplace_back(f); }, num_threads);
    std::string input_decoded;
    no_init_resize(input_decoded, input.size());
    lz77_sss::decode(factorization.begin(), input_decoded.data(), input.size());
    EXPECT_EQ(input, input_decoded);

    if (exact) {
        std::vector<uint64_t> lengths;
        for (const lz77_sss::factor& f : factorization) lengths.push_back(f.text_len());
        EXPECT_EQ(lengths, reference_phrase_lengths(input)) << "n=" << input.size() << " threads=" << num_threads;
    }
}

std::string random_fasta_input(uint64_t max_size) {
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    const bool full_alphabet = prob(gen) < 0.5;
    const std::string src = random_repetitive_input<std::string>(1, max_size,
        full_alphabet ? std::numeric_limits<char>::min() : 'A',
        full_alphabet ? std::numeric_limits<char>::max() : 'Z');
    const double header_rate = prob(gen);
    std::uniform_int_distribution<uint64_t> line_len(0, random_log_uniform_size(1, 1000, gen));
    std::string input;
    uint64_t i = 0;

    while (i < src.size()) {
        const uint64_t len = std::min<uint64_t>(src.size() - i, line_len(gen));
        if (prob(gen) < header_rate) input.push_back(prob(gen) < 0.5 ? '>' : ';');
        input.append(src, i, len);
        i += len;
        if (i < src.size() || prob(gen) < 0.5) input.push_back('\n');
    }

    return input;
}

void test_fasta_roundtrip(const std::string& input) {
    uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    const text_encoding encoding = text_encoding(std::uniform_int_distribution<int>(0, 3)(gen));
    const bool exact = std::uniform_int_distribution<int>(0, 1)(gen) == 1;
    const fasta_mode mode = fasta_mode(std::uniform_int_distribution<int>(0, 2)(gen));
    const std::string path = (std::filesystem::temp_directory_path() /
        ("lz77_sss_fasta_" + random_alphanumeric_string(10))).string();
    std::ofstream(path, std::ios::binary).write(input.data(), input.size());
    lz77_sss::parameters params { .num_threads = num_threads,
        .transf_mode = std::uniform_int_distribution<int>(0, 1)(gen) == 1 ? lz77_sss::with_interval_samples : lz77_sss::without_interval_samples };
    std::vector<lz77_sss::factor> factorization;
    auto sink = [&](lz77_sss::factor f) { factorization.emplace_back(f); };
    fasta_headers headers;

    auto factorize = [&](auto text, auto out) {
        if (exact) lz77_sss::factorize_exact(text, out, params);
        else lz77_sss::factorize_approximate(text, out, params);
    };

    with_text_from_file(path, input.size(), encoding, exact ? exact_factorization : aprx_factorization, mode, headers,
        4 * lz77_sss::default_tau, num_threads, false, [&](auto T) {
        fasta_interleaver interleaver(std::move(headers), sink, [&](char* text, uint64_t size, auto out) {
            factorize(lz77_sss::direct_text(text, size), out);
        });

        factorize(T, [&](lz77_sss::factor f) { interleaver.add(f); });
        interleaver.finish();
    });

    std::filesystem::remove(path);
    std::string input_decoded;
    no_init_resize(input_decoded, input.size());
    lz77_sss::decode(factorization.begin(), input_decoded.data(), input.size());
    EXPECT_EQ(input, input_decoded);
}

void test_file_decoder(const std::string& input) {
    const uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    const bool in_ram = std::uniform_int_distribution<int>(0, 1)(gen) == 1;
    const uint64_t buffer_bits = std::uniform_int_distribution<uint64_t>(2, 20)(gen);
    std::string text = input;
    std::vector<lz77_sss::factor> factorization;
    run_approximate(text, [&](lz77_sss::factor f) { factorization.emplace_back(f); }, num_threads);
    const std::string path = (std::filesystem::temp_directory_path() /
        ("lz77_sss_decode_" + random_alphanumeric_string(10))).string();

    {
        file_decoder decoder(path, input.size(), in_ram, num_threads, buffer_bits);

        for (const lz77_sss::factor& f : factorization) {
            const char c = char(uint64_t(f.src));

            if (!f.is_literal()) {
                decoder.copy(f.src, f.len);
            } else if (std::uniform_int_distribution<int>(0, 1)(gen) == 0) {
                decoder.literal(c);
            } else {
                decoder.literals(&c, 1);
            }
        }

        decoder.finish();
        EXPECT_TRUE(decoder.good());
    }

    std::string output;
    no_init_resize(output, input.size());
    std::ifstream(path, std::ios::binary).read(output.data(), input.size());
    EXPECT_EQ(std::filesystem::file_size(path), input.size());
    std::filesystem::remove(path);
    EXPECT_EQ(input, output) << "n=" << input.size() << " in_ram=" << in_ram << " buffer_bits=" << buffer_bits;
}

template <typename int_t>
struct int_input {
    std::vector<uint32_t> ranks;
    std::vector<uint64_t> values;
    std::vector<int_t> input;
};

template <typename int_t>
int_input<int_t> random_int_input(uint64_t max_size) {
    const uint64_t n = random_log_uniform_size(1, max_size, gen);
    uint64_t max_bits = 40;
    if constexpr (!std::is_same_v<int_t, uint40_t>) max_bits = 8 * sizeof(int_t);
    const uint64_t sigma = random_log_uniform_size(1, std::clamp<uint64_t>(n / 4, 1,
        uint64_t { 1 } << std::min<uint64_t>(20, max_bits)), gen);
    int_input<int_t> result;
    std::vector<uint32_t>& ranks = result.ranks;
    ranks = random_repetitive_input<std::vector<uint32_t>>(n, n, uint32_t(0), uint32_t(sigma - 1));
    std::uniform_real_distribution<double> prob(0.0, 1.0);

    if (prob(gen) < 0.5) {
        const double frequent_share = prob(gen);
        std::uniform_int_distribution<uint32_t> frequent(0, uint32_t(std::min<uint64_t>(sigma, random_log_uniform_size(1, 64, gen)) - 1));

        for (uint32_t& c : ranks) {
            if (prob(gen) < frequent_share) c = frequent(gen);
        }
    }

    constexpr uint64_t mersenne_prime = (uint64_t { 1 } << 61) - 1;

    if (max_bits == 64 && prob(gen) < 0.25) {
        for (uint64_t c = 0; c < sigma; c++) result.values.push_back((c % 8) * mersenne_prime + c / 8);
    } else {
        const uint64_t range_bits = std::uniform_int_distribution<uint64_t>(std::bit_width(sigma - 1), max_bits)(gen);
        const uint64_t value_range = range_bits == 64 ? UINT64_MAX : std::max<uint64_t>(sigma, uint64_t { 1 } << range_bits);
        std::unordered_set<uint64_t> chosen;

        for (uint64_t j = value_range - sigma; j < value_range; j++) {
            uint64_t v = std::uniform_int_distribution<uint64_t>(0, j)(gen);

            if (!chosen.insert(v).second) {
                v = j;
                chosen.insert(j);
            }

            result.values.push_back(v);
        }
    }

    std::shuffle(result.values.begin(), result.values.end(), gen);
    result.input.resize(n);
    for (uint64_t i = 0; i < n; i++) result.input[i] = int_t(result.values[ranks[i]]);
    return result;
}

template <typename int_t>
std::vector<uint64_t> to_u64(const std::vector<int_t>& v) {
    return std::vector<uint64_t>(v.begin(), v.end());
}

std::vector<uint64_t> reference_phrase_lengths(const std::vector<uint32_t>& input) {
    std::vector<uint64_t> lengths;
    lz77::emit_function emit = [&](lz77::factor f) { lengths.push_back(f.text_len()); };
    lz77::parallel_lpf_factorizer().factorize(input.begin(), input.end(), emit, emit, 1);
    return lengths;
}

template <typename int_t>
void test_int_roundtrip_of(uint64_t max_size, lz77_sss::parameters params, bool exact) {
    int_input<int_t> inp = random_int_input<int_t>(max_size);
    std::vector<int_t>& input = inp.input;
    const std::vector<uint64_t> original = to_u64(input);
    const uint64_t n = input.size();
    params.num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    std::vector<lz77_sss::factor> factorization;
    auto out = [&](lz77_sss::factor f) { factorization.emplace_back(f); };
    const int kind = std::uniform_int_distribution<int>(0, 3)(gen);

    if (kind == 0) {
        if (exact) lz77_sss::factorize_exact(input.data(), n, out, params);
        else lz77_sss::factorize_approximate(input.data(), n, out, params);
        EXPECT_EQ(to_u64(input), original) << "n=" << n << " bytes=" << sizeof(int_t);
    } else {
        auto run = [&](auto T) {
            auto mapped = [&](lz77_sss::factor f) {
                if (f.is_literal()) f.src = inp.values[f.src];
                out(f);
            };

            if (exact) lz77_sss::factorize_exact(T, mapped, params);
            else lz77_sss::factorize_approximate(T, mapped, params);
        };

        const uint64_t sigma = inp.values.size();
        if (kind == 1) run(lz77_sss::int_direct_text(inp.ranks.data(), n, sigma));
        else if (kind == 2) run(lz77_sss::int_packed_text(inp.ranks.data(), n, sigma, params.num_threads));
        else run(lz77_sss::int_split_text(inp.ranks.data(), n, sigma, params.num_threads));
    }

    std::vector<int_t> decoded(n);
    lz77_sss::decode(factorization.begin(), decoded.data(), n);
    EXPECT_EQ(to_u64(decoded), original) << "n=" << n << " kind=" << kind << " bytes=" << sizeof(int_t);

    if (exact) {
        std::vector<uint64_t> lengths;
        for (const lz77_sss::factor& f : factorization) lengths.push_back(f.text_len());
        EXPECT_EQ(lengths, reference_phrase_lengths(inp.ranks)) << "n=" << n << " kind=" << kind
            << " bytes=" << sizeof(int_t) << " threads=" << params.num_threads;
    }
}

void test_int_roundtrip(uint64_t max_size, lz77_sss::parameters params, bool exact) {
    switch (std::uniform_int_distribution<int>(0, 8)(gen)) {
        case 0: test_int_roundtrip_of<int8_t>(max_size, params, exact); break;
        case 1: test_int_roundtrip_of<uint8_t>(max_size, params, exact); break;
        case 2: test_int_roundtrip_of<int16_t>(max_size, params, exact); break;
        case 3: test_int_roundtrip_of<uint16_t>(max_size, params, exact); break;
        case 4: test_int_roundtrip_of<int32_t>(max_size, params, exact); break;
        case 5: test_int_roundtrip_of<uint32_t>(max_size, params, exact); break;
        case 6: test_int_roundtrip_of<uint40_t>(max_size, params, exact); break;
        case 7: test_int_roundtrip_of<int64_t>(max_size, params, exact); break;
        default: test_int_roundtrip_of<uint64_t>(max_size, params, exact); break;
    }
}

lz77_sss::parameters random_exact_params(lz77_sss::transform_mode transf_mode) {
    const range_ds_type type = std::uniform_int_distribution<int>(0, 1)(gen) == 0 ? range_ds_type::sdsg : range_ds_type::swkdt;
    return { .fact_mode = lz77_sss::auto_gaps, .transf_mode = transf_mode,
        .range_ds = { .type = type, .decomposed = true }, .exact_alg = lz77_sss::sss_based };
}

TEST(test_lz77_sss, approximate) {
    run_fuzz("lz77-sss", {
        { "approximate", [](uint64_t) { test_roundtrip(random_repetitive_input<std::string>(1, 500000), run_approximate); }, false },
    }, fuzz_iterations(1500));
}

TEST(test_lz77_sss, approximate_exact_gaps) {
    run_fuzz("lz77-sss", {
        { "approximate-exact-gaps", [](uint64_t) { test_roundtrip(random_repetitive_input<std::string>(1, 500000), run_approximate_exact_gaps); }, false },
    }, fuzz_iterations(1500));
}

TEST(test_lz77_sss, exact_with_interval_samples) {
    run_fuzz("lz77-sss", {
        { "exact-with-interval-samples", [](uint64_t) { test_roundtrip(random_repetitive_input<std::string>(1, 100000), run_exact_with_interval_samples, true); }, false },
    }, fuzz_iterations(1000));
}

TEST(test_lz77_sss, exact_without_interval_samples) {
    run_fuzz("lz77-sss", {
        { "exact-without-interval-samples", [](uint64_t) { test_roundtrip(random_repetitive_input<std::string>(1, 100000), run_exact_without_interval_samples, true); }, false },
    }, fuzz_iterations(1000));
}

TEST(test_lz77_sss, exact_sa) {
    run_fuzz("lz77-sss", {
        { "exact-sa", [](uint64_t) { test_roundtrip(random_repetitive_input<std::string>(1, 100000), run_exact_sa, true); }, false },
    }, fuzz_iterations(1500));
}

TEST(test_lz77_sss, exact_auto) {
    run_fuzz("lz77-sss", {
        { "exact-auto", [](uint64_t) { test_roundtrip(random_repetitive_input<std::string>(1, 100000), run_exact_auto, true); }, false },
    }, fuzz_iterations(1000));
}

TEST(test_lz77_sss, file_decoder) {
    run_fuzz("lz77-sss", {
        { "file-decoder", [](uint64_t) { test_file_decoder(random_repetitive_input<std::string>(1, 200000)); }, false },
    }, fuzz_iterations(1000));
}

void test_fasta_reader(const std::string& input) {
    const uint16_t num_threads = std::uniform_int_distribution<uint16_t>(1, omp_get_max_threads())(gen);
    const uint64_t block = random_log_uniform_size(16, 4096, gen);
    const std::string path = (std::filesystem::temp_directory_path() /
        ("lz77_sss_fasta_reader_" + random_alphanumeric_string(10))).string();
    std::ofstream(path, std::ios::binary).write(input.data(), input.size());

    fasta_headers expected_headers;
    std::string expected;
    no_init_resize(expected, input.size());
    fasta_filter filter(&expected_headers);
    expected.resize(filter.filter(input.data(), input.size(), expected.data(), 0));
    if (expected_headers.text_end.size() < expected_headers.seq_pos.size()) expected_headers.text_end.push_back(expected_headers.text.size());

    fasta_headers headers;
    const std::optional<fasta_layout> layout = scan_fasta(path, input.size(), headers,
        std::numeric_limits<uint64_t>::max(), num_threads, [](const char*, uint64_t, uint16_t) { }, block);
    ASSERT_TRUE(layout.has_value());
    std::string sequence;
    no_init_resize(sequence, layout->seq_off.back());
    read_fasta_sequence(path, input.size(), *layout, num_threads,
        [&](const char* data, uint64_t len, uint64_t at, uint16_t) { std::memcpy(sequence.data() + at, data, len); });
    EXPECT_EQ(sequence, expected) << "n=" << input.size() << " block=" << block << " threads=" << num_threads;
    EXPECT_EQ(headers.text, expected_headers.text);
    EXPECT_EQ(headers.seq_pos, expected_headers.seq_pos);
    EXPECT_EQ(headers.text_end, expected_headers.text_end);
    EXPECT_EQ(headers.bytes_removed, input.size() - expected.size());

    if (!expected.empty()) {
        char_histogram histogram { };
        for (char c : expected) histogram[uint8_t(c)]++;
        lce::text::packed_text<> text(histogram, expected.size());
        parallel_packer<lce::text::packed_text<>> packer(text, num_threads);
        read_fasta_sequence(path, input.size(), *layout, num_threads,
            [&](const char* data, uint64_t len, uint64_t at, uint16_t t) { packer.pack(data, len, at, t); });
        packer.finish();
        uint64_t mismatches = 0;
        for (uint64_t i = 0; i < expected.size(); i++) mismatches += text.char_at(i) != expected[i];
        EXPECT_EQ(mismatches, uint64_t { 0 }) << "n=" << input.size() << " block=" << block << " threads=" << num_threads;
    }

    std::filesystem::remove(path);
}

TEST(test_lz77_sss, fasta_reader) {
    run_fuzz("lz77-sss", {
        { "fasta-reader", [](uint64_t) { test_fasta_reader(random_fasta_input(100000)); }, false },
    }, fuzz_iterations(1400));
}

TEST(test_lz77_sss, fasta) {
    run_fuzz("lz77-sss", {
        { "fasta", [](uint64_t) { test_fasta_roundtrip(random_fasta_input(50000)); }, false },
    }, fuzz_iterations(1000));
}

template <typename int_t>
void test_int_extreme_values(std::vector<int_t> input) {
    const std::vector<uint64_t> original = to_u64(input);

    for (bool exact : { false, true }) {
        std::vector<lz77_sss::factor> factorization;
        auto out = [&](lz77_sss::factor f) { factorization.emplace_back(f); };
        if (exact) lz77_sss::factorize_exact(input.data(), input.size(), out);
        else lz77_sss::factorize_approximate(input.data(), input.size(), out);
        std::vector<int_t> decoded(input.size());
        lz77_sss::decode(factorization.begin(), decoded.data(), decoded.size());
        EXPECT_EQ(to_u64(decoded), original) << "exact=" << exact;
        EXPECT_EQ(to_u64(input), original) << "exact=" << exact;
    }
}

TEST(test_lz77_sss, int_extreme_values) {
    constexpr uint64_t mersenne_prime = (uint64_t { 1 } << 61) - 1;
    test_int_extreme_values<uint64_t>({ UINT64_MAX, 0, mersenne_prime, UINT64_MAX, 0, mersenne_prime, 7, UINT64_MAX });
    test_int_extreme_values<int64_t>({ -1, 7, INT64_MIN, -1, 7, INT64_MAX, INT64_MIN, -1, 7 });
    test_int_extreme_values<int32_t>({ -1, INT32_MIN, INT32_MAX, -1, INT32_MIN, 0 });
}

TEST(test_lz77_sss, approximate_int) {
    run_fuzz("lz77-sss", {
        { "approximate-int", [](uint64_t) { test_int_roundtrip(200000, { .fact_mode = lz77_sss::auto_gaps }, false); }, false },
    }, fuzz_iterations(1500));
}

TEST(test_lz77_sss, approximate_exact_gaps_int) {
    run_fuzz("lz77-sss", {
        { "approximate-exact-gaps-int", [](uint64_t) { test_int_roundtrip(200000, { .fact_mode = lz77_sss::exact_gaps }, false); }, false },
    }, fuzz_iterations(1500));
}

TEST(test_lz77_sss, exact_with_interval_samples_int) {
    run_fuzz("lz77-sss", {
        { "exact-with-interval-samples-int", [](uint64_t) { test_int_roundtrip(50000, random_exact_params(lz77_sss::with_interval_samples), true); }, false },
    }, fuzz_iterations(1000));
}

TEST(test_lz77_sss, exact_without_interval_samples_int) {
    run_fuzz("lz77-sss", {
        { "exact-without-interval-samples-int", [](uint64_t) { test_int_roundtrip(50000, random_exact_params(lz77_sss::without_interval_samples), true); }, false },
    }, fuzz_iterations(1000));
}

TEST(test_lz77_sss, exact_sa_int) {
    run_fuzz("lz77-sss", {
        { "exact-sa-int", [](uint64_t) { test_int_roundtrip(50000, { .exact_alg = lz77_sss::sa_based }, true); }, false },
    }, fuzz_iterations(1500));
}

TEST(test_lz77_sss, exact_auto_int) {
    run_fuzz("lz77-sss", {
        { "exact-auto-int", [](uint64_t) { test_int_roundtrip(50000, { .exact_alg = lz77_sss::auto_select }, true); }, false },
    }, fuzz_iterations(1000));
}
