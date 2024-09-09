#include <fstream>

std::string text_name;
std::string result_file_path;
#define LZ77_SSS_BENCH 1

#include <lz77_sss/lz77_sss.hpp>

int main(int argc, char** argv) {
    if (!(2 <= argc && argc <= 5)) {
        std::cout << "usage: lz77_sss_bench-tau <file> <result_file> <min_tau> <max_tau>";
        std::cout << "       the last three parameters are optional";
        std::cout << "       min_tau and max_tau must be in the range [4, 4096]";
        exit(-1);
    }

    std::string file_path = argv[1];
    std::ifstream input_file(file_path);
    
    if (!input_file.good()) {
        std::cout << "error: could not read <file>";
        exit(-1);
    }

    if (argc >= 3) {
        text_name = file_path.substr(
            file_path.find_last_of("/\\") + 1);
        result_file_path = argv[2];
    }

    uint64_t min_tau = 4;
    uint64_t max_tau = 4096;

    if (argc >= 4) min_tau = atol(argv[3]);
    if (argc >= 5) max_tau = atol(argv[4]);

    input_file.seekg(0, std::ios::end);
    uint64_t n = input_file.tellg();
    input_file.seekg(0, std::ios::beg);
    auto t0 = now();
    std::cout << "reading T (" << format_size(n) << ")" << std::flush;
    std::string T;
    no_init_resize(T, n);
    read_from_file(input_file, T.data(), n);
    input_file.close();
    log_runtime(t0);

    for_constexpr_pow<4, 4096>([&](auto tau){
        if (min_tau <= tau && tau <= max_tau) {
            std::cout << std::endl <<
                "running LZ77 SSS 3-approximation with tau = " <<
                tau << ":" << std::endl;
            std::ofstream fact_sss_file("fact_sss_aprx");
            
            if (n <= std::numeric_limits<uint32_t>::max()) {
                lz77_sss<uint32_t>::factorize_approximate<
                    greedy, lpf_all_external, tau>(T, fact_sss_file, true);
            } else {
                lz77_sss<uint64_t>::factorize_approximate<
                    greedy, lpf_all_external, tau>(T, fact_sss_file, true);
            }

            fact_sss_file.close();
            std::filesystem::remove("fact_sss_aprx");
        }
    });
}