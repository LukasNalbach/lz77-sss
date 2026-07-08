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

#include <lz77_sss/lz77_sss.hpp>

int main(int argc, char** argv)
{
    if (!(2 <= argc && argc <= 5)) {
        std::cout << "usage: lz77_sss_bench-tau <file> <result_file> <min_tau> <max_tau>" << std::endl;
        std::cout << "       the last three parameters are optional" << std::endl;
        std::cout << "       min_tau and max_tau must be in the range [4, 4096]" << std::endl;
        exit(-1);
    }

    std::string file_path = argv[1];
    std::ifstream input_file(file_path);

    if (!input_file.good()) {
        std::cout << "error: could not read <file>" << std::endl;
        exit(-1);
    }

    if (argc >= 3) {
        text_name = file_path.substr(
            file_path.find_last_of("/\\") + 1);
        result_file_path = argv[2];
        result_file.open(result_file_path, std::ofstream::app);
    }

    uint64_t min_tau = 4;
    uint64_t max_tau = 4096;

    if (argc >= 4)
        min_tau = atol(argv[3]);
    if (argc >= 5)
        max_tau = atol(argv[4]);

    uint64_t n = std::filesystem::file_size(argv[1]);
    auto t0 = now();
    std::cout << "reading input (" << format_size(n) << ")" << std::flush;
    std::string T;
    no_init_resize_with_excess(T, n, 4 * lz77_sss<>::default_tau);
    read_from_file(input_file, T.data(), n);
    input_file.close();
    log_runtime(t0);

    for_constexpr_pow<4, 4096>([&](auto tau) {
        if (min_tau <= tau && tau <= max_tau) {
            std::cout << std::endl <<
                "running LZ77 SSS 3-approximation with tau = "
                << tau << ":" << std::endl;
            std::ofstream fact_sss_file("fact_sss_aprx");

            if (n <= std::numeric_limits<uint32_t>::max()) {
                lz77_sss<uint32_t>::factorize_approximate<
                    greedy, lpf_opt, tau>(T.data(), n,
                        [&](auto f){fact_sss_file << f;},
                        { .num_threads = 1, .log = true });
            } else {
                lz77_sss<uint64_t>::factorize_approximate<
                    greedy, lpf_opt, tau>(T.data(), n,
                        [&](auto f){fact_sss_file << f;},
                        { .num_threads = 1, .log = true });
            }

            fact_sss_file.close();
            std::filesystem::remove("fact_sss_aprx");
        }
    });
    
    return 0;
}