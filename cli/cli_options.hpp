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

#include <filesystem>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <optional>
#include <string>

#include <lz77_sss/misc/text_reader.hpp>

struct cli_options {
    std::string input_file_path;
    std::string output_file_path;
    uint16_t num_threads = 1;
    text_encoding encoding = auto_encoding;
    fasta_mode fasta = fasta_auto;
    range_ds_kind range_ds = lz77_sss::parameters { }.range_ds;
    lz77_sss::exact_algorithm exact_alg = lz77_sss::parameters { }.exact_alg;
};

inline void help(std::string name, bool exact, std::string message)
{
    if (message != "") std::cout << message << std::endl;
    std::cout << "usage: " << name << " [options] <input_file>" << std::endl;
    std::cout << " -o <output_file>  output file (default: <input_file>" << (exact ? ".lz77-exact" : ".lz77-aprx")
              << ")" << std::endl;
    std::cout << " -t <threads>      number of threads (default: all)" << std::endl;
    std::cout << " -enc <encoding>   how the text is kept in memory (default: auto):" << std::endl;
    std::cout << "                   plain   one byte per character (fastest)" << std::endl;
    std::cout << "                   packed  fewer bits per character for small alphabets" << std::endl;
    std::cout << "                   split   fewer bits for frequent characters (slowest)" << std::endl;
    std::cout << "                   auto    packed or split if that saves memory, else plain" << std::endl;
    std::cout << " -fasta <mode>     handle the header lines of FASTA files separately from" << std::endl;
    std::cout << "                   the sequences: on, off or auto (default: auto, on for" << std::endl;
    std::cout << "                   FASTA files)" << std::endl;

    if (exact) {
        std::cout << " -alg <algorithm>  exact algorithm (default: auto):" << std::endl;
        std::cout << "                   sss     string synchronizing sets, for repetitive texts" << std::endl;
        std::cout << "                   sa      suffix array, needs about 9 bytes per character" << std::endl;
        std::cout << "                   auto    sa if it needs at most "
                  << std::round(100 * (lz77_sss::max_sa_peak_ratio - 1)) << " % more memory than sss," << std::endl;
        std::cout << "                           else sss" << std::endl;
        std::cout << " -r <range_ds>     range data structure (default: d-swsg):" << std::endl;
        std::cout << "                   swsg    static weighted square grid" << std::endl;
        std::cout << "                   swkdt   static weighted k-d tree" << std::endl;
        std::cout << "                   sdsg    semi-dynamic square grid (single-threaded)" << std::endl;
        std::cout << "                   d-...   one structure per character, e.g. d-swsg" << std::endl;
    }

    std::cout << " -h                show help" << std::endl;
    exit(-1);
}

inline cli_options parse_cli_args(int argc, char** argv, std::string name, bool exact = false)
{
    cli_options options;
    options.num_threads = omp_get_max_threads();
    if (argc == 1) help(name, exact, "");
    int arg_idx = 1;

    while (arg_idx < argc - 1) {
        std::string arg = argv[arg_idx++];

        if (arg == "-o") {
            if (arg_idx >= argc - 1) help(name, exact, "error: missing parameter after -o option");
            options.output_file_path = argv[arg_idx++];
        } else if (arg == "-t") {
            if (arg_idx >= argc - 1) help(name, exact, "error: missing parameter after -t option");
            options.num_threads = std::max<uint16_t>(1, atoi(argv[arg_idx++]));
            if (options.num_threads > omp_get_max_threads()) help(name, exact, "error: requested too many threads");
        } else if (arg == "-enc") {
            if (arg_idx >= argc - 1) help(name, exact, "error: missing parameter after -enc option");
            options.encoding = parse_text_encoding(argv[arg_idx++], [&](std::string message) { help(name, exact, message); });
        } else if (arg == "-fasta") {
            if (arg_idx >= argc - 1) help(name, exact, "error: missing parameter after -fasta option");
            options.fasta = parse_fasta_mode(argv[arg_idx++], [&](std::string message) { help(name, exact, message); });
        } else if (exact && arg == "-alg") {
            if (arg_idx >= argc - 1) help(name, exact, "error: missing parameter after -alg option");
            std::string alg = argv[arg_idx++];

            if (alg == "auto") {
                options.exact_alg = lz77_sss::auto_select;
            } else if (alg == "sss") {
                options.exact_alg = lz77_sss::sss_based;
            } else if (alg == "sa") {
                options.exact_alg = lz77_sss::sa_based;
            } else {
                help(name, exact, "error: unknown exact algorithm '" + alg + "'");
            }
        } else if (exact && arg == "-r") {
            if (arg_idx >= argc - 1) help(name, exact, "error: missing parameter after -r option");
            std::string range_ds = argv[arg_idx++];
            std::optional<range_ds_kind> kind = range_ds_by_name(range_ds);
            if (!kind) help(name, exact, "error: unknown range data structure '" + range_ds + "'");
            options.range_ds = *kind;
        } else if (arg == "-h") {
            help(name, exact, "");
        } else {
            help(name, exact, "error: unrecognized '" + arg + "' option");
        }
    }

    options.input_file_path = argv[argc - 1];
    if (options.input_file_path == "-h") help(name, exact, "");
    if (options.input_file_path[0] == '-') help(name, exact, "error: missing <input_file>");
    if (options.output_file_path == "") options.output_file_path = options.input_file_path + (exact ? ".lz77-exact" : ".lz77-aprx");
    return options;
}
