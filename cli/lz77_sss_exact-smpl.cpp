#include <fstream>
#include <lz77_sss/lz77_sss.hpp>

int main(int argc, char** argv)
{
    if (!(3 <= argc && argc <= 4)) {
        std::cout << "usage: lz77_sss_exact-smpl <input_file> <output_file> <threads>" << std::endl;
        std::cout << "       the last parameter is optional" << std::endl;
        exit(-1);
    }

    std::ifstream input_file(argv[1]);
    std::ofstream output_file(argv[2]);
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
    output_file << n;
    std::cout << std::endl << "running LZ77 SSS exact algorithm (with samples):" << std::endl;

    if (n <= std::numeric_limits<uint32_t>::max()) {
        lz77_sss<uint32_t>::factorize_exact<
            greedy, lpf_opt, with_samples>(T,
            output_file, { .num_threads = p, .log = true });
    } else {
        lz77_sss<uint64_t>::factorize_exact<
            greedy, lpf_opt, with_samples>(T,
            output_file, { .num_threads = p, .log = true });
    }

    output_file.close();
}