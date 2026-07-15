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

std::string text_name;
std::string result_file_path;
std::ofstream result_file;
#define LZ77_SSS_BENCH 1

#include <lz77/lpf.hpp>
#include <lz77/kkp2.hpp>
#include <lz77_sss/lz77_sss.hpp>

uint16_t max_threads = 0;
std::string T;
uint64_t n;

template <
    factorize_mode fact_mode,
    phrase_mode phr_mode>
void run_sss_approximate(std::string file_name, uint16_t max_threads)
{

    for (uint16_t num_threads = 1; num_threads <= max_threads; num_threads *= 2) {
        std::ofstream fact_sss_file(file_name);

        if (n <= std::numeric_limits<uint32_t>::max()) {
            lz77_sss<uint32_t>::factorize_approximate<
                fact_mode, phr_mode>(T.data(), n,
                    [&](auto f){fact_sss_file << f;},
                    { .num_threads = num_threads, .log = true });
        } else {
            lz77_sss<uint64_t>::factorize_approximate<
                fact_mode, phr_mode>(T.data(), n,
                    [&](auto f){fact_sss_file << f;},
                    { .num_threads = num_threads, .log = true });
        }
    }
}

template <
    factorize_mode fact_mode,
    phrase_mode phr_mode,
    transform_mode transf_mode,
    template <typename> typename range_ds_t>
void run_sss_exact(std::string file_name, uint16_t max_threads)
{

    for (uint16_t num_threads = 1; num_threads <= max_threads; num_threads *= 2) {
        std::ofstream fact_sss_file(file_name);

        if (n <= std::numeric_limits<uint32_t>::max()) {
            lz77_sss<uint32_t>::factorize_exact<
                fact_mode, phr_mode, transf_mode, range_ds_t>(
                T.data(), n, [&](auto f){fact_sss_file << f;},
                { .num_threads = num_threads, .log = true });
        } else {
            lz77_sss<uint64_t>::factorize_exact<
                fact_mode, phr_mode, transf_mode, range_ds_t>(
                T.data(), n, [&](auto f){fact_sss_file << f;},
                { .num_threads = num_threads, .log = true });
        }
    }
}

void log_algorithm(
    std::string alg_name, uint64_t time,
    uint64_t mem_peak, uint64_t num_factors,
    uint16_t num_threads = 1)
{
    double comp_ratio = n / (double)num_factors;

    std::cout << ", in ~ " << format_time(time) << std::endl;
    std::cout << "throughput: " << format_throughput(n, time) << std::endl;
    std::cout << "peak memory consumption: " << format_size(mem_peak) << std::endl;
    std::cout << "input length / num. of factors: " << comp_ratio << std::endl;

    if (result_file_path != "") {
        result_file << "RESULT"
            << " text_name=" << text_name
            << " n=" << n
            << " alg=" << alg_name
            << " num_threads=" << num_threads
            << " num_factors=" << num_factors
            << " comp_ratio=" << comp_ratio
            << " time=" << time
            << " throughput=" << throughput(n, time)
            << " mem_peak=" << mem_peak << std::endl;
    }
}

int main(int argc, char** argv)
{
    if (!(2 <= argc && argc <= 4)) {
        std::cout << "usage: lz77_sss_bench <file> <max_threads> <result_file>" << std::endl;
        std::cout << "       the last two parameters are optional" << std::endl;
        exit(-1);
    }

    std::string file_path = argv[1];
    std::ifstream input_file(file_path);

    if (!input_file.good()) {
        std::cout << "error: could not read <file>" << std::endl;
        exit(-1);
    }

    if (argc >= 3) {
        max_threads = atoi(argv[2]);

        if (max_threads == 0 || max_threads > omp_get_max_threads()) {
            std::cout << "error: invalid number of threads" << std::endl;
            exit(-1);
        }
    }

    if (argc == 4) {
        text_name = file_path.substr(
            file_path.find_last_of("/\\") + 1);
        result_file_path = argv[3];
        result_file.open(result_file_path, std::ofstream::app);
    }

    n = std::filesystem::file_size(argv[1]);
    auto t0 = now();
    std::cout << "reading input (" << format_size(n) << ")" << std::flush;
    no_init_resize_with_excess(T, n, 4 * lz77_sss<>::default_tau);
    read_from_file(input_file, T.data(), n);
    input_file.close();
    log_runtime(t0);

    std::cout << std::endl << "running LZ77 SSS 3-approximation:" << std::endl;
    run_sss_approximate<greedy, lpf_opt>("fact_sss_aprx", max_threads);
    std::filesystem::remove("fact_sss_aprx");
    
    std::cout << std::endl << "running naive LZ77 SSS exact algorithm:" << std::endl;
    run_sss_exact<greedy, lpf_opt, naive,
        static_weighted_square_grid>("fact_sss_exact", 1);
    std::filesystem::remove("fact_sss_exact");

    std::cout << std::endl << "running LZ77 SSS exact algorithm (without samples):" << std::endl;
    run_sss_exact<greedy, lpf_opt, without_samples,
        decomposed_static_weighted_square_grid>("fact_sss_exact", max_threads);
    std::filesystem::remove("fact_sss_exact");
    
    std::cout << std::endl << "running LZ77 SSS exact algorithm (with samples):" << std::endl;
    run_sss_exact<greedy, lpf_opt, with_samples,
        decomposed_static_weighted_square_grid>("fact_sss_exact", max_threads);
    std::filesystem::remove("fact_sss_exact");

    for (uint16_t num_threads = 1; num_threads <= max_threads; num_threads *= 2) {
        std::cout << std::endl << "running LZ77 LPF algorithm (" << num_threads << " threads)" << std::flush;
        uint64_t baseline_memory_alloc = malloc_count_current();
        std::ofstream file_lpf("fact_lpf");
        malloc_count_reset_peak();
        auto t1 = now();
        lz77::parallel_lpf_factorizer().factorize(T.begin(), T.end(),
            std::ostream_iterator<lz77::factor>(file_lpf, ""), num_threads, "fact_lpf_tmp");
        auto t2 = now();
        uint64_t num_factors = file_lpf.tellp() / sizeof(lz77::factor);
        uint64_t time = time_diff_ns(t1, t2);
        uint64_t mem_peak = malloc_count_peak() - baseline_memory_alloc;
        file_lpf.close();
        std::filesystem::remove("fact_lpf");
        log_algorithm("lpf", time, mem_peak, num_factors, num_threads);
    }

    for (uint16_t num_threads = 1; num_threads <= max_threads; num_threads *= 2) {
        std::cout << std::endl << "running LZ77 KKP2 algorithm (" << num_threads << " threads)" << std::flush;
        uint64_t baseline_memory_alloc = malloc_count_current();
        std::ofstream file_kkp2("fact_kkp2");
        malloc_count_reset_peak();
        auto t1 = now();
        lz77::parallel_kkp2_factorizer().factorize(T.begin(), T.end(),
            std::ostream_iterator<lz77::factor>(file_kkp2, ""), num_threads, "fact_kkp2_tmp");
        auto t2 = now();
        uint64_t num_factors = file_kkp2.tellp() / sizeof(lz77::factor);
        uint64_t time = time_diff_ns(t1, t2);
        uint64_t mem_peak = malloc_count_peak() - baseline_memory_alloc;
        file_kkp2.close();
        std::filesystem::remove("fact_kkp2");
        log_algorithm("kkp2", time, mem_peak, num_factors, num_threads);
    }

    return 0;
}