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
#include <lz77_sss/misc/file_decoder.hpp>
#include <lz77_sss/misc/huffman.hpp>

void help(std::string message)
{
    if (message != "") std::cout << message << std::endl;
    std::cout << "usage: lz77-sss-decode [options] <input_file> <output_file>" << std::endl;
    std::cout << " -ram          keep the whole output in memory while decoding (default: only" << std::endl;
    std::cout << "               the last 64 MiB, older parts are read back from <output_file>)" << std::endl;
    std::cout << " -t <threads>  number of threads (default: all)" << std::endl;
    std::cout << " -h            show help" << std::endl;
    exit(-1);
}

int main(int argc, char** argv)
{
    bool in_ram = false;
    uint16_t num_threads = omp_get_max_threads();
    int arg_idx = 1;

    while (arg_idx < argc - 2) {
        std::string arg = argv[arg_idx++];

        if (arg == "-ram") {
            in_ram = true;
        } else if (arg == "-t") {
            if (arg_idx >= argc - 2) help("error: missing parameter after -t option");
            num_threads = std::max<uint16_t>(1, atoi(argv[arg_idx++]));
            if (num_threads > omp_get_max_threads()) help("error: requested too many threads");
        } else if (arg == "-h") {
            help("");
        } else {
            help("error: unrecognized '" + arg + "' option");
        }
    }

    if (argc - arg_idx != 2) help("");
    const std::string input_file_path = argv[arg_idx];
    const std::string output_file_path = argv[arg_idx + 1];
    std::ifstream input_file(input_file_path, std::ios::in | std::ios::binary);

    if (!input_file.good()) {
        std::cout << "error: could not read <input_file>" << std::endl;
        exit(-1);
    }

    uint64_t input_file_size = std::filesystem::file_size(input_file_path);
    std::cout << "input file size = " << format_size(input_file_size) << std::endl;
    uint64_t n = 0;
    input_file.read((char*) &n, 5);
    file_decoder decoder(output_file_path, n, in_ram, num_threads);

    if (!decoder.good()) {
        std::cout << "error: could not write to <output_file>" << std::endl;
        exit(-1);
    }

    bit_reader reader(input_file);
    huffman len_huff, dist_huff;
    using factor = lz77_sss::factor;
    huff_factor_iterator<factor> it(reader, len_huff, dist_huff, n);
    log_phase_begin(true, "decoding (" + format_size(n) + ", " + (in_ram ? "in main memory" : "low-space") + ")");
    auto t1 = now();

    {
        phase_progress progress(true, n);
        uint64_t next_report = progress.grain();

        while (decoder.position() < n) {
            const factor f = *it++;

            if (f.is_literal()) {
                decoder.literal(char(uint64_t(f.src)));
            } else {
                decoder.copy(f.src, f.len);
            }

            if (decoder.position() >= next_report) [[unlikely]] {
                progress.reached(decoder.position());
                next_report = decoder.position() + progress.grain();
            }
        }

        decoder.finish();
    }

    auto t2 = now();
    log_runtime(t1, t2);

    if (!decoder.good()) {
        std::cout << "error: writing <output_file> failed" << std::endl;
        exit(-1);
    }

    std::cout << "throughput = " << format_throughput(n, time_diff_ns(t1, t2)) << std::endl;
    std::cout << "peak memory consumption = " << format_size(malloc_count_peak()) << std::endl;
    std::cout << "compression ratio = " << n / (double) input_file_size << std::endl;
    return 0;
}
