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
#include <filesystem>
#include <ips4o.hpp>

#include <lz77_sss/data_structures/range/range.hpp>
#include <lz77_sss/misc/log.hpp>

struct point {
    uint64_t x;
    uint64_t y;
    uint64_t weight;
};

struct __attribute__((packed)) operation {
    bool is_insert;
    char c;
    uint64_t weight;
    uint64_t x_1;
    uint64_t x_2;
    uint64_t y_1;
    uint64_t y_2;
};

std::string T;
uint64_t n;
std::vector<uint64_t> sampling;
std::vector<point> points;
std::vector<operation> operations;

uint64_t min_win_size = 1 << 11;
uint64_t max_win_size = 1 << 16;

void help(std::string message)
{
    if (message != "") std::cout << message << std::endl;
    std::cout << "usage: bench-range-queries [-m <m_file>] <input_file> <queries_file>" << std::endl;
    std::cout << "                           [<min_win> <max_win>]" << std::endl;
    std::cout << " -m <m_file>     append results to <m_file>" << std::endl;
    std::cout << " <queries_file>  queries written by gen-range-queries for <input_file>" << std::endl;
    std::cout << " <min_win>       log2 of the smallest square grid window (default: 11)" << std::endl;
    std::cout << " <max_win>       log2 of the largest square grid window (default: 16)" << std::endl;
    exit(-1);
}

void bench(const range_ds_kind& kind, uint64_t win_size) {
    std::cout << "benchmarking " << kind.name() << std::flush;
    if (win_size != 0) std::cout << " with window size " << win_size << std::flush;

    bench_win_size = win_size;
    using point_t = range_ds::point_t;
    std::vector<uint40_t> sampling_local;
    std::vector<point_t> points_local;

    for (uint64_t s : sampling) {
        sampling_local.emplace_back(s);
    }

    uint64_t num_iterations = 0;
    uint64_t time_build = 0;
    uint64_t time_query = 0;
    uint64_t mem_peak = 0;
    uint64_t mem_used = 0;
    uint64_t check_sum;
    auto time_start = now();

    while (time_diff_sec(time_start, now()) < 2) {
        num_iterations++;
        points_local.clear();
        points_local.reserve(points.size());

        for (point p : points) {
            point_t p_add {};
            p_add.x = p.x;
            p_add.y = p.y;
            if (kind.is_static()) p_add.weight = p.weight;
            points_local.emplace_back(p_add);
        }

        range_ds* ds = nullptr;
        malloc_count_reset_peak();
        uint64_t baseline_bytes = malloc_count_current();
        auto t0 = now();

        if (kind.decomposed) {
            ds = make_range_ds(kind, T.data(), sampling_local, points_local, 1);
        } else {
            ds = make_range_ds(kind, points_local, points.size(), 1);
        }

        auto t1 = now();

        if (kind.is_static()) {
            points_local.clear();
            points_local.shrink_to_fit();
        }

        check_sum = 0;
        time_build += time_diff_ns(t0, t1);
        t1 = now();

        for (operation op : operations) {
            if (op.is_insert) {
                if (kind.is_dynamic()) {
                    point_t p {};
                    p.x = op.x_1;
                    p.y = op.y_1;
                    ds->insert(op.c, p);
                }
            } else {
                point_t p;
                bool result;

                if (kind.is_static()) {
                    std::tie(p, result) = ds->lighter_point_in_range(
                        op.c, op.weight,
                        op.x_1, op.x_2,
                        op.y_1, op.y_2);
                } else {
                    std::tie(p, result) = ds->point_in_range(
                        op.c,
                        op.x_1, op.x_2,
                        op.y_1, op.y_2);
                }

                if (result) {
                    check_sum += result;
                    check_sum += op.x_1;
                    check_sum += op.x_2;
                    check_sum += op.y_1;
                    check_sum += op.y_2;
                    check_sum += op.weight;
                }
            }
        }

        auto t2 = now();
        time_query += time_diff_ns(t1, t2);
        mem_used = ds->size_in_bytes();
        mem_peak = malloc_count_peak() < baseline_bytes ?
            0 : (malloc_count_peak() - baseline_bytes);
        delete ds;
    }

    time_query /= num_iterations;
    time_build /= num_iterations;
    uint64_t num_points = points.size();
    uint64_t num_operations = operations.size();
    uint64_t num_queries = std::count_if(operations.begin(), operations.end(),
        [](const operation& op) { return !op.is_insert; });

    std::cout << std::endl;
    std::cout << "construction time = " << time_build / (1.0 * num_points) << " ns/point" << std::endl;
    std::cout << "construction memory peak = " << mem_peak / (1.0 * num_points) << " bytes/point" << std::endl;
    std::cout << "throughput = " << (num_operations * 1000.0) / time_query << " inserts & queries/us" << std::endl;
    std::cout << "size = " << mem_used / (1.0 * num_points) << " bytes/point" << std::endl;
    std::cout << "checksum = " << check_sum << std::endl;
    std::cout << std::endl;

    if (result_log::out.is_open()) {
        result_log::out << "RESULT"
            << " text_name=" << result_log::text_name
            << " ds=" << kind.name()
            << " win_size=" << win_size
            << " num_points=" << num_points
            << " num_queries=" << num_queries
            << " num_operations=" << num_operations
            << " time_build=" << time_build
            << " time_query=" << time_query
            << " mem_peak=" << mem_peak
            << " mem_used=" << mem_used
            << std::endl;
    }
}

void bench_all() {
    for (bool decomposed : { false, true }) {
        for (uint64_t win_size = min_win_size; win_size <= max_win_size; win_size *= 2) {
            bench(range_ds_kind { .type = range_ds_type::sdsg, .decomposed = decomposed }, win_size);
        }
    }

    for (bool decomposed : { false, true }) {
        for (uint64_t win_size = min_win_size; win_size <= max_win_size; win_size *= 2) {
            bench(range_ds_kind { .type = range_ds_type::swsg, .decomposed = decomposed }, win_size);
        }
    }

    bench(range_ds_kind { .type = range_ds_type::swkdt }, 0);
    bench(range_ds_kind { .type = range_ds_type::swkdt, .decomposed = true }, 0);
}

int main(int argc, char** argv)
{
    int arg_idx = 1;

    if (arg_idx + 1 < argc && std::string(argv[arg_idx]) == "-m") {
        result_log::path = argv[arg_idx + 1];
        arg_idx += 2;
    }

    if ((argc - arg_idx != 2 && argc - arg_idx != 4) || argv[arg_idx][0] == '-') help("");
    omp_set_num_threads(1);
    std::string file_path = argv[arg_idx];
    result_log::text_name = file_path.substr(file_path.find_last_of("/\\") + 1);
    std::ifstream input_file(file_path);
    std::ifstream queries_file(argv[arg_idx + 1]);
    if (!input_file.good()) help("error: could not read <input_file>");
    if (!queries_file.good()) help("error: could not read <queries_file>");

    if (argc - arg_idx == 4) {
        const uint64_t log2_min_win_size = atoi(argv[arg_idx + 2]);
        const uint64_t log2_max_win_size = atoi(argv[arg_idx + 3]);
        if (log2_min_win_size > log2_max_win_size || log2_max_win_size > 32) help("error: invalid range of window sizes");
        min_win_size = uint64_t { 1 } << log2_min_win_size;
        max_win_size = uint64_t { 1 } << log2_max_win_size;
    }

    if (result_log::path != "") {
        result_log::out.open(result_log::path, std::ios::app);
        if (!result_log::out.good()) help("error: could not write to <m_file>");
    }

    n = std::filesystem::file_size(file_path);
    auto t0 = now();
    std::cout << "reading input (" << format_size(n) << ")" << std::flush;
    no_init_resize_with_excess(T, n, 4 * 4096);
    read_fully(input_file, T.data(), n);
    input_file.close();
    log_runtime(t0);

    uint64_t num_points;
    queries_file.read((char*) &num_points, 8);
    sampling.reserve(num_points);
    points.reserve(num_points);

    std::cout << "reading samples, points and operations" << std::flush;

    for (uint64_t i = 0; i < num_points; i++) {
        uint64_t s;
        queries_file.read((char*) &s, 8);
        sampling.emplace_back(s);
    }

    for (uint64_t i = 0; i < num_points; i++) {
        point p;
        queries_file.read((char*) &p, sizeof(point));
        points.emplace_back(p);
    }

    while (queries_file.peek() != EOF) {
        operation op;
        queries_file.read((char*) &op, sizeof(operation));
        operations.emplace_back(op);
    }

    std::cout << std::endl << std::endl;

    bench_all();
    return 0;
}
