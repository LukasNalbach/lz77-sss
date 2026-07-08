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
#include <lz77_sss/lz77_sss.hpp>
#include <lz77_sss/misc/huffman.hpp>

int main(int argc, char** argv)
{
    if (!(3 <= argc && argc <= 4)) {
        std::cout << "usage: lz77_sss_exact <input_file> <output_file> <threads>" << std::endl;
        std::cout << "       the last parameter is optional" << std::endl;
        exit(-1);
    }

    std::ifstream input_file(argv[1], std::ios::binary);
    std::ofstream output_file(argv[2], std::ios::binary);
    uint16_t p = omp_get_max_threads();

    if (!input_file.good()) {
        std::cout << "error: could not read <input_file>" << std::endl;
        exit(-1);
    }

    if (!output_file.good()) {
        std::cout << "error: could not write to <output_file>" << std::endl;
        exit(-1);
    }

    if (argc >= 4) {
        p = atoi(argv[3]);

        if (p == 0 || p > omp_get_max_threads()) {
            std::cout << "error: invalid number of threads" << std::endl;
            exit(-1);
        }
    }

    uint64_t n = std::filesystem::file_size(argv[1]);
    auto t0 = now();
    std::cout << "reading input (" << format_size(n) << ")" << std::flush;
    std::string T;
    no_init_resize_with_excess(T, n, 4 * lz77_sss<>::default_tau);
    read_from_file(input_file, T.data(), n);
    input_file.close();
    log_runtime(t0);
    std::cout << "running LZ77 SSS exact algorithm (with samples):" << std::endl;
    huff_writer writer(output_file, n);

    if (n <= std::numeric_limits<uint32_t>::max()) {
        lz77_sss<uint32_t>::factorize_exact<
            greedy, lpf_opt, with_samples>(T.data(), n,
            [&](auto f) { writer.add(f); },
            { .num_threads = p, .log = true });
    } else {
        lz77_sss<uint64_t>::factorize_exact<
            greedy, lpf_opt, with_samples>(T.data(), n,
            [&](auto f) { writer.add(f); },
            { .num_threads = p, .log = true });
    }

    writer.finish();
    output_file.close();
    uint64_t output_file_size = std::filesystem::file_size(argv[2]);
    std::cout << "output file size: " << format_size(output_file_size) << std::endl;
    std::cout << "compression ratio: " << n / (double) output_file_size << std::endl;
    return 0;
}