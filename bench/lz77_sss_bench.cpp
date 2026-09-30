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

#include <fstream>

#include <lz77/parallel_lpf.hpp>
#include <lz77/kkp2.hpp>
#include <lz77_sss/inst.hpp>
#include <lz77_sss/misc/text_reader.hpp>

uint64_t n;

void help(std::string message)
{
    if (message != "") std::cout << message << std::endl;
    std::cout << "usage: lz77-sss-bench [-m <m_file>] <input_file> [<max_threads>]" << std::endl;
    std::cout << " -m <m_file>    append results to <m_file>" << std::endl;
    std::cout << " <max_threads>  largest thread count; the count doubles from 1 (default: all)" << std::endl;
    exit(-1);
}

template <typename text_t>
void run_sss_approximate(const text_t& T, lz77_sss::factorize_mode fact_mode,
    std::string file_name, uint16_t max_threads)
{
    for (uint16_t num_threads = 1; num_threads <= max_threads; num_threads *= 2) {
        std::ofstream fact_sss_file(file_name);

        lz77_sss::factorize_approximate(T,
                [&](auto f){fact_sss_file << f;},
                { .num_threads = num_threads, .log = true,
                  .fact_mode = fact_mode });
    }
}

template <typename text_t>
void run_sss_exact(const text_t& T, lz77_sss::factorize_mode fact_mode,
    lz77_sss::transform_mode transf_mode, const range_ds_kind& range_ds,
    std::string file_name, uint16_t max_threads)
{
    for (uint16_t num_threads = 1; num_threads <= max_threads; num_threads *= 2) {
        std::ofstream fact_sss_file(file_name);

        lz77_sss::factorize_exact(
            T, [&](auto f){fact_sss_file << f;},
            { .num_threads = num_threads, .log = true,
              .fact_mode = fact_mode, .transf_mode = transf_mode,
              .range_ds = range_ds, .exact_alg = lz77_sss::sss_based });
    }
}

void log_algorithm(
    std::string alg_name, uint64_t time,
    uint64_t mem_peak, uint64_t num_factors,
    uint16_t num_threads = 1)
{
    double comp_ratio = n / (double)num_factors;

    std::cout << ", in ~ " << format_time(time) << std::endl;
    std::cout << "throughput = " << format_throughput(n, time) << std::endl;
    std::cout << "peak memory consumption = " << format_size(mem_peak) << std::endl;
    std::cout << "input length / num. of factors = " << comp_ratio << std::endl;

    if (result_log::out.is_open()) {
        result_log::out << "RESULT"
            << " text_name=" << result_log::text_name
            << " n=" << n
            << " alg=" << alg_name
            << " num_threads=" << num_threads
            << " num_factors=" << num_factors
            << " comp_ratio=" << comp_ratio
            << " time=" << time
            << " throughput=" << throughput_mb_per_s(n, time)
            << " mem_peak=" << mem_peak << std::endl;
    }
}

int main(int argc, char** argv)
{
    int arg_idx = 1;

    if (arg_idx + 1 < argc && std::string(argv[arg_idx]) == "-m") {
        result_log::path = argv[arg_idx + 1];
        arg_idx += 2;
    }

    if (!(1 <= argc - arg_idx && argc - arg_idx <= 2) || argv[arg_idx][0] == '-') help("");
    std::string file_path = argv[arg_idx];
    std::ifstream input_file(file_path);
    if (!input_file.good()) help("error: could not read <input_file>");
    uint16_t max_threads = omp_get_max_threads();

    if (argc - arg_idx == 2) {
        max_threads = atoi(argv[arg_idx + 1]);
        if (max_threads == 0 || max_threads > omp_get_max_threads()) help("error: invalid number of threads");
    }

    if (result_log::path != "") {
        if (!std::ofstream(result_log::path, std::ofstream::app).good()) help("error: could not write to <m_file>");
        result_log::text_name = file_path.substr(file_path.find_last_of("/\\") + 1);
    }

    n = std::filesystem::file_size(file_path);
    input_file.close();
    fasta_headers headers;

    auto open_result_log = []() {
        if (result_log::path == "") return;
        result_log::out.open(result_log::path, std::ofstream::app);
        result_log::write_rows = true;
    };

    with_text_from_file(file_path, n, auto_encoding, aprx_factorization, fasta_off, headers,
        4 * lz77_sss::default_tau, omp_get_max_threads(), true, [&](auto T) {
        open_result_log();
        std::cout << std::endl << "running LZ77 SSS 3-approximation:" << std::endl;
        run_sss_approximate(T, lz77_sss::auto_gaps, "fact_sss_aprx", max_threads);
        std::filesystem::remove("fact_sss_aprx");
    });

    if (result_log::out.is_open()) result_log::out.close();
    std::cout << std::endl;

    with_text_from_file(file_path, n, auto_encoding, exact_factorization, fasta_off, headers,
        4 * lz77_sss::default_tau, omp_get_max_threads(), true, [&](auto T) {
        open_result_log();
        std::cout << std::endl << "running LZ77 SSS exact algorithm (without interval samples):" << std::endl;
        run_sss_exact(T, lz77_sss::auto_gaps, lz77_sss::without_interval_samples,
            { .type = range_ds_type::swsg, .decomposed = true }, "fact_sss_exact", max_threads);
        std::filesystem::remove("fact_sss_exact");

        std::cout << std::endl << "running LZ77 SSS exact algorithm (with interval samples):" << std::endl;
        run_sss_exact(T, lz77_sss::auto_gaps, lz77_sss::with_interval_samples,
            { .type = range_ds_type::swsg, .decomposed = true }, "fact_sss_exact", max_threads);
        std::filesystem::remove("fact_sss_exact");
    });

    std::string T;
    auto t0 = now();
    std::cout << std::endl << "reading input for the reference algorithms ("
        << format_size(n) << ")" << std::flush;
    no_init_resize_with_excess(T, n, 4 * lz77_sss::default_tau);
    direct_ifstream reference_file(file_path);
    read_fully(reference_file, T.data(), n);
    reference_file.close();
    log_runtime(t0);

    for (uint16_t num_threads = 1; num_threads <= max_threads; num_threads *= 2) {
        std::cout << std::endl << "running LZ77 LPF algorithm (" << num_threads << " threads)" << std::flush;
        uint64_t baseline_bytes = malloc_count_current();
        std::ofstream file_lpf("fact_lpf");
        malloc_count_reset_peak();
        auto t1 = now();
        lz77::parallel_lpf_factorizer().factorize(T.begin(), T.end(),
            std::ostream_iterator<lz77::factor>(file_lpf, ""), num_threads, "fact_lpf_tmp");
        auto t2 = now();
        uint64_t num_factors = file_lpf.tellp() / sizeof(lz77::factor);
        uint64_t time = time_diff_ns(t1, t2);
        uint64_t mem_peak = malloc_count_peak() - baseline_bytes;
        file_lpf.close();
        std::filesystem::remove("fact_lpf");
        log_algorithm("lpf", time, mem_peak, num_factors, num_threads);
    }

    std::cout << std::endl << "running LZ77 KKP2 algorithm" << std::flush;
    uint64_t baseline_bytes = malloc_count_current();
    std::ofstream file_kkp2("fact_kkp2");
    malloc_count_reset_peak();
    auto t1 = now();
    lz77::kkp2_factorizer().factorize(T.begin(), T.end(),
        std::ostream_iterator<lz77::factor>(file_kkp2, ""));
    auto t2 = now();
    uint64_t num_factors = file_kkp2.tellp() / sizeof(lz77::factor);
    uint64_t time = time_diff_ns(t1, t2);
    uint64_t mem_peak = malloc_count_peak() - baseline_bytes;
    file_kkp2.close();
    std::filesystem::remove("fact_kkp2");
    log_algorithm("kkp2", time, mem_peak, num_factors);

    return 0;
}
