#include <fstream>
#include <ips4o.hpp>

uint64_t cur_win_size;
#define BENCH_RANGE_QUERIES 1

#include <lz77_sss/data_structures/static_weighted_range/static_weighted_kd_tree.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_square_grid.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_striped_square.hpp>

#include <lz77_sss/data_structures/dynamic_range/dynamic_square_grid.hpp>
#include <lz77_sss/data_structures/dynamic_range/semi_dynamic_square_grid.hpp>

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

const char* T = nullptr;
uint64_t n;
std::vector<uint64_t> sampling;
std::vector<point> points;
std::vector<operation> operations;
std::vector<operation> queries;

std::ofstream results_file;
std::string text_name;
uint64_t min_win_size = 1 << 11;
uint64_t max_win_size = 1 << 16;

template <typename pos_t, typename sidx_t, template <typename> typename range_ds_t>
void bench(uint64_t win_size, std::string log_name) {
    std::cout << "benchmarking " <<  range_ds_t<sidx_t>::name() << std::flush;
    if (win_size != 0) std::cout << " with window size " << win_size << std::flush;

    cur_win_size = win_size;
    using point_t = range_ds_t<sidx_t>::point_t;
    std::vector<pos_t> sampling_local;
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
            point_t p_add { .x = p.x, .y = p.y };
            if constexpr (range_ds_t<pos_t>::is_static()) p_add.weight = p.weight;
            points_local.emplace_back(p_add);
        }

        range_ds_t<sidx_t> ds;
        malloc_count_reset_peak();
        uint64_t alloc_before = malloc_count_current();
        auto t0 = now();

        if constexpr (range_ds_t<sidx_t>::is_decomposed()) {
            ds = range_ds_t<sidx_t>(T, sampling_local, points_local);
        } else {
            ds = range_ds_t<sidx_t>(points_local, points.size());
        }

        auto t1 = now();
        
        if constexpr (range_ds_t<sidx_t>::is_static()) {
            points_local.clear();
            points_local.shrink_to_fit();
        }
        
        check_sum = 0;
        time_build += time_diff_ns(t0, t1);
        t1 = now();

        for (operation op : operations) {
            if (op.is_insert) {
                if constexpr (range_ds_t<sidx_t>::is_dynamic()) {
                    point_t p { .x = op.x_1, .y = op.y_1 };

                    if constexpr (range_ds_t<sidx_t>::is_decomposed()) {
                        ds.insert(op.c, p);
                    } else {
                        ds.insert(p);
                    }
                }
            } else {
                point_t p;
                bool result;

                if constexpr (range_ds_t<sidx_t>::is_static()) {
                    if constexpr (range_ds_t<sidx_t>::is_decomposed()) {
                        std::tie(p, result) = ds.lighter_point_in_range(
                            op.c, op.weight,
                            op.x_1, op.x_2,
                            op.y_1, op.y_2);
                    } else {
                        std::tie(p, result) = ds.lighter_point_in_range(
                            op.weight,
                            op.x_1, op.x_2,
                            op.y_1, op.y_2);
                    }
                } else {
                    if constexpr (range_ds_t<sidx_t>::is_decomposed()) {
                        std::tie(p, result) = ds.point_in_range(
                            op.c,
                            op.x_1, op.x_2,
                            op.y_1, op.y_2);
                    } else {
                        std::tie(p, result) = ds.point_in_range(
                            op.x_1, op.x_2,
                            op.y_1, op.y_2);
                    }
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
        mem_used = ds.size_in_bytes();
        mem_peak = malloc_count_peak() < alloc_before ?
            0 : (malloc_count_peak() - alloc_before);
    }

    time_query /= num_iterations;
    time_build /= num_iterations;
    uint64_t num_points = points.size();
    uint64_t num_operations = operations.size();
    uint64_t num_queries = num_operations - num_points;

    std::cout << std::endl;
    std::cout << "construction time: " << time_build / (1.0 * num_points) << " ns/point" << std::endl;
    std::cout << "construction memory peak: " << mem_peak / (1.0 * num_points) << " bytes/point" << std::endl;
    std::cout << (num_operations * 1000.0) / time_query << " inserts & queries/us" << std::endl;
    std::cout << "size: " << mem_used / (1.0 * num_points) << " bytes/point" << std::endl;
    std::cout << "checksum: " << check_sum << std::endl;
    std::cout << std::endl;

    if (results_file.is_open()) {
        results_file << "RESULT"
            << " text="  << text_name
            << " ds=" << log_name
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

template <typename pos_t, typename sidx_t>
void bench_all() {
    for (uint64_t win_size = min_win_size; win_size <= max_win_size; win_size *= 2) {
        bench<pos_t, sidx_t, dynamic_square_grid>(win_size, "dsg");
    }

    for (uint64_t win_size = min_win_size; win_size <= max_win_size; win_size *= 2) {
        bench<pos_t, sidx_t, decomposed_dynamic_square_grid>(win_size, "ddsg");
    }

    for (uint64_t win_size = min_win_size; win_size <= max_win_size; win_size *= 2) {
        bench<pos_t, sidx_t, semi_dynamic_square_grid>(win_size, "sdsg");
    }

    for (uint64_t win_size = min_win_size; win_size <= max_win_size; win_size *= 2) {
        bench<pos_t, sidx_t, decomposed_semi_dynamic_square_grid>(win_size, "dsdsg");
    }

    for (uint64_t win_size = min_win_size; win_size <= max_win_size; win_size *= 2) {
        bench<pos_t, sidx_t, static_weighted_square_grid>(win_size, "swsg");
    }

    for (uint64_t win_size = min_win_size; win_size <= max_win_size; win_size *= 2) {
        bench<pos_t, sidx_t, decomposed_static_weighted_square_grid>(win_size, "dswsg");
    }
    
    bench<pos_t, sidx_t, static_weighted_striped_square>(0, "swss");
    bench<pos_t, sidx_t, decomposed_static_weighted_striped_square>(0, "dswss");
    bench<pos_t, sidx_t, static_weighted_kd_tree>(0, "swkdt");
    bench<pos_t, sidx_t, decomposed_static_weighted_kd_tree>(0, "dswkdt");
}

int main(int argc, char** argv)
{
    if (!(3 <= argc && argc <= 6)) {
        std::cout << "usage: bench_range_queries <file> <queries_file> <results_file> <min_win_size> <max_win_size>" << std::endl;
        std::cout << "       the last three parameters are optional" << std::endl;
        exit(-1);
    }

    omp_set_num_threads(1);
    std::string file_path = argv[1];
    text_name = file_path.substr(file_path.find_last_of("/\\") + 1);
    std::ifstream input_file(file_path);
    std::ifstream queries_file(argv[2]);

    if (!input_file.good()) {
        std::cout << "error: could not read <input_file>" << std::endl;
        exit(-1);
    }

    if (!queries_file.good()) {
        std::cout << "error: could not read <queries_file>" << std::endl;
        exit(-1);
    }

    if (argc >= 4) {
        results_file.open(argv[3], std::ios::app);

        if (!results_file.good()) {
            std::cout << "error: could not read <results_file>" << std::endl;
            exit(-1);
        }
    }

    if (argc >= 5) min_win_size = 1 << atoi(argv[4]);
    if (argc >= 6) max_win_size = 1 << atoi(argv[5]);
    
    input_file.seekg(0, std::ios::end);
    n = input_file.tellg();
    input_file.seekg(0, std::ios::beg);
    auto t0 = now();
    std::cout << "reading T (" << format_size(n) << ")" << std::flush;
    char* T = (char*) std::aligned_alloc(16, n + 4 * 4096);
    read_from_file(input_file, T, n);
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

    if (n <= std::numeric_limits<uint32_t>::max()) {
        bench_all<uint32_t, uint32_t>();
    } else {
        if (num_points <= std::numeric_limits<uint32_t>::max()) {
            bench_all<uint64_t, uint32_t>();
        } else {
            bench_all<uint64_t, uint64_t>();
        }
    }
}