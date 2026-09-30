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

#include <lz77_sss/inst.hpp>
#include <lz77_sss/misc/text_reader.hpp>

void help(std::string message)
{
    if (message != "") std::cout << message << std::endl;
    std::cout << "usage: lz77-sss-bench-tau [-m <m_file>] <input_file> [<min_tau> <max_tau>]" << std::endl;
    std::cout << " -m <m_file>  append results to <m_file>" << std::endl;
    std::cout << " <min_tau>    smallest tau (default: 4, at least 4)" << std::endl;
    std::cout << " <max_tau>    largest tau (default: 4096, at most 4096); every power of two" << std::endl;
    std::cout << "              in between is run" << std::endl;
    exit(-1);
}

int main(int argc, char** argv)
{
    int arg_idx = 1;

    if (arg_idx + 1 < argc && std::string(argv[arg_idx]) == "-m") {
        result_log::path = argv[arg_idx + 1];
        arg_idx += 2;
    }

    if ((argc - arg_idx != 1 && argc - arg_idx != 3) || argv[arg_idx][0] == '-') help("");
    std::string file_path = argv[arg_idx];
    std::ifstream input_file(file_path);
    if (!input_file.good()) help("error: could not read <input_file>");
    uint64_t min_tau = 4;
    uint64_t max_tau = 4096;

    if (argc - arg_idx == 3) {
        min_tau = atol(argv[arg_idx + 1]);
        max_tau = atol(argv[arg_idx + 2]);
    }

    if (min_tau < 4 || max_tau > 4096 || min_tau > max_tau) help("error: invalid range of tau");

    if (result_log::path != "") {
        if (!std::ofstream(result_log::path, std::ofstream::app).good()) help("error: could not write to <m_file>");
        result_log::text_name = file_path.substr(file_path.find_last_of("/\\") + 1);
    }

    uint64_t n = std::filesystem::file_size(file_path);
    input_file.close();
    fasta_headers headers;

    with_text_from_file(file_path, n, auto_encoding, aprx_factorization, fasta_off, headers,
        4 * max_tau, omp_get_max_threads(), true, [&](auto T) {
        if (result_log::path != "") {
            result_log::out.open(result_log::path, std::ofstream::app);
            result_log::write_rows = true;
        }

        for (uint64_t tau = std::bit_ceil(min_tau); tau <= max_tau; tau *= 2) {
            std::cout << std::endl <<
                "running LZ77 SSS 3-approximation with tau = "
                << tau << ":" << std::endl;
            std::ofstream fact_sss_file("fact_sss_aprx");

            lz77_sss::factorize_approximate(T,
                    [&](auto f){fact_sss_file << f;},
                    { .num_threads = 1, .log = true, .tau = tau,
                      .fact_mode = lz77_sss::auto_gaps });

            fact_sss_file.close();
            std::filesystem::remove("fact_sss_aprx");
        }
    });

    return 0;
}
