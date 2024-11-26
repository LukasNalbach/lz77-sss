#include <fstream>
#include <lz77_sss/lz77_sss.hpp>

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cout << "usage: lz77_sss_decode <input_file> <output_file>" << std::endl;
        exit(-1);
    }

    std::ifstream input_file(argv[1]);
    std::ofstream output_file(argv[2]);

    if (!input_file.good()) {
        std::cout << "error: could not read <input_file>" << std::endl;
        exit(-1);
    }

    if (!output_file.good()) {
        std::cout << "error: could not write to <output_file>" << std::endl;
        exit(-1);
    }

    uint64_t n;
    input_file >> n;
    std::cout << "decoding T (" << format_size(n) << ")" << std::flush;
    std::string T;
    auto time = now();

    if (n <= std::numeric_limits<uint32_t>::max()) {
        T = lz77_sss<uint32_t>::decode(input_file, n);
    } else {
        T = lz77_sss<uint64_t>::decode(input_file, n);
    }

    input_file.close();
    time = log_runtime(time);
    std::cout << "writing T to the output file" << std::flush;
    write_to_file(output_file, T.data(), n);
    time = log_runtime(time);
    output_file.close();
}