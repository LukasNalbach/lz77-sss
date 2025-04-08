#include <fstream>

std::ofstream queries_file;
#define GEN_RANGE_QUERIES 1
#include <lz77_sss/lz77_sss.hpp>

int main(int argc, char** argv)
{
    if (!(3 <= argc && argc <= 4)) {
        std::cout << "usage: gen_range_queries <file> <queries_file>" << std::endl;
        exit(-1);
    }

    std::ifstream input_file(argv[1]);
    queries_file.open(argv[2]);

    if (!input_file.good()) {
        std::cout << "error: could not read <file>" << std::endl;
        exit(-1);
    }

    if (!queries_file.good()) {
        std::cout << "error: could not read <queries_file>" << std::endl;
        exit(-1);
    }

    input_file.seekg(0, std::ios::end);
    uint64_t n = input_file.tellg();
    input_file.seekg(0, std::ios::beg);
    auto t0 = now();
    std::cout << "reading T (" << format_size(n) << ")" << std::flush;
    std::string T;
    no_init_resize_with_excess(T, n, 4 * lz77_sss<>::default_tau);
    read_from_file(input_file, T.data(), n);
    input_file.close();
    log_runtime(t0);

    std::cout << "generating queries" << std::flush;
    std::ofstream fact_sss_file("fact_sss_exact");

    if (n <= std::numeric_limits<uint32_t>::max()) {
        lz77_sss<uint32_t>::factorize_exact<
            greedy, lpf_opt, without_samples, decomposed_static_weighted_kd_tree>(
            T.data(), n, [&](auto f){fact_sss_file << f;}, { .num_threads = 1, .log = false });
    } else {
        lz77_sss<uint64_t>::factorize_exact<
            greedy, lpf_opt, without_samples, decomposed_static_weighted_kd_tree>(
            T.data(), n, [&](auto f){fact_sss_file << f;}, { .num_threads = 1, .log = false });
    }

    std::filesystem::remove("fact_sss_exact");
    std::cout << std::endl;
}