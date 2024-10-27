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

    input_file.seekg(0, std::ios::end);
    uint64_t n = input_file.tellg();
    input_file.seekg(0, std::ios::beg);
    auto t0 = now();
    std::cout << "reading T (" << format_size(n) << ")" << std::flush;
    std::string T;
    no_init_resize_with_exess(T, n, 4 * 4096);
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
                    greedy, lpf_opt, tau>(T,
                    fact_sss_file, { .num_threads = 1, .log = true });
            } else {
                lz77_sss<uint64_t>::factorize_approximate<
                    greedy, lpf_opt, tau>(T,
                    fact_sss_file, { .num_threads = 1, .log = true });
            }

            fact_sss_file.close();
            std::filesystem::remove("fact_sss_aprx");
        }
    });
}