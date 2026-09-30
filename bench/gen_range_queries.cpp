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

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cout << "usage: gen-range-queries <input_file> <queries_file>" << std::endl;
        exit(-1);
    }

    std::ifstream input_file(argv[1]);
    result_log::queries_out.open(argv[2], std::ios::binary);

    if (!input_file.good()) {
        std::cout << "error: could not read <input_file>" << std::endl;
        exit(-1);
    }

    if (!result_log::queries_out.good()) {
        std::cout << "error: could not write to <queries_file>" << std::endl;
        exit(-1);
    }

    uint64_t n = std::filesystem::file_size(argv[1]);
    input_file.close();
    fasta_headers headers;

    with_text_from_file(argv[1], n, auto_encoding, exact_factorization, fasta_off, headers,
        4 * lz77_sss::default_tau, omp_get_max_threads(), true, [&](auto T) {
        std::cout << "generating queries" << std::flush;
        std::ofstream fact_sss_file("fact_sss_exact");

        lz77_sss::factorize_exact(
            T, [&](auto f){fact_sss_file << f;},
            { .num_threads = 1, .log = false,
              .fact_mode = lz77_sss::auto_gaps, .transf_mode = lz77_sss::without_interval_samples,
              .range_ds = { .type = range_ds_type::swkdt, .decomposed = true },
              .exact_alg = lz77_sss::sss_based });
    });

    result_log::queries_out.close();
    std::filesystem::remove("fact_sss_exact");
    std::cout << std::endl;
    return 0;
}
