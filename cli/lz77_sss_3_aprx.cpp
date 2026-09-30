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
#include <lz77_sss/misc/huffman.hpp>
#include <lz77_sss/misc/fasta_interleaver.hpp>

#include "cli_options.hpp"

int main(int argc, char** argv)
{
    auto time_start = now();
    cli_options options = parse_cli_args(argc, argv, "lz77-sss-3-aprx");
    direct_ifstream input_file(options.input_file_path);
    direct_ofstream output_file(options.output_file_path);

    if (!input_file.good()) {
        std::cout << "error: could not read " << options.input_file_path << std::endl;
        exit(-1);
    }

    if (!output_file.good()) {
        std::cout << "error: could not write to " << options.output_file_path << std::endl;
        exit(-1);
    }

    uint64_t n = std::filesystem::file_size(options.input_file_path);
    fasta_headers headers;
    input_file.close();

    with_text_from_file(options.input_file_path, n, options.encoding, aprx_factorization, options.fasta, headers,
        4 * lz77_sss::default_tau, options.num_threads, true, [&](auto T) {
        std::cout << "running LZ77 SSS 3-approximation:" << std::endl;
        huff_factor_writer writer(output_file, n);
        auto sink = [&](auto f) { writer.add(f); };
        lz77_sss::parameters params { .num_threads = options.num_threads, .fact_mode = lz77_sss::auto_gaps };

        fasta_interleaver interleaver(std::move(headers), sink, [&](char* text, uint64_t size, auto out) {
            lz77_sss::factorize_approximate(text, size, out, params);
        }, true);

        params.log = true;
        params.time_start = time_start;
        params.log_factor_count = false;
        lz77_sss::factorize_approximate(T, [&](auto f) { interleaver.add(f); }, params);
        interleaver.finish();
        writer.finish();

        std::cout << "num. of factors = " << interleaver.num_factors() << std::endl;
        std::cout << "input length / num. of factors = "
                  << n / (double) std::max<uint64_t>(1, interleaver.num_factors()) << std::endl;
    });

    output_file.close();

    uint64_t output_file_size = std::filesystem::file_size(options.output_file_path);
    std::cout << "output file size = " << format_size(output_file_size) << std::endl;
    std::cout << "compression ratio = " << n / (double) output_file_size << std::endl;
    return 0;
}