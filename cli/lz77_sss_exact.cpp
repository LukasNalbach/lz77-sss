#include <fstream>
#include <lz77_sss/lz77_sss.hpp>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cout << "usage: lz77_sss_exact <input_file> <output_file>";
        exit(-1);
    }

    std::ifstream input_file(argv[1]);
    std::ofstream output_file(argv[2]);
    
    if (!input_file.good()) {
        std::cout << "error: could not read <input_file>";
        exit(-1);
    }
    
    if (!output_file.good()) {
        std::cout << "error: could not write to <output_file>";
        exit(-1);
    }
    
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
    output_file << n;
    std::cout << std::endl << "running LZ77 SSS exact algorithm (without samples):" << std::endl;

    if (n <= std::numeric_limits<uint32_t>::max()) {
        lz77_sss<uint32_t>::factorize_exact<
            greedy, lpf_all, without_samples>(T, output_file, true);
    } else {
        lz77_sss<uint64_t>::factorize_approximate<
            greedy, lpf_all, without_samples>(T, output_file, true);
    }

    output_file.close();
}