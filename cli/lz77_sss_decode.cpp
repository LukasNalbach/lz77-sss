#include <fstream>
#include <lz77_sss/lz77_sss.hpp>

uint64_t n;
std::fstream input_file;
std::fstream output_file;

template <typename pos_t>
void decode()
{
    using factor = lz77_sss<pos_t>::factor;
    factor f;
    pos_t pos_input = 5;
    pos_t pos_output = 0;
    uint64_t buff_size = std::max<uint64_t>(32 * 1024, n / 1000);
    std::string buff;

    while (pos_output < n) {
        input_file >> f;
        pos_input += factor::size_of();

        if (f.len == 0) {
            output_file.write((char*) &f.src, 1);
            pos_input++;
            pos_output++;
        } else {
            copy_buffered(output_file, output_file,
                buff, f.src, pos_output, f.len, buff_size);
            pos_output += f.len;
        }
    }
}

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cout << "usage: lz77_sss_decode <input_file> <output_file>" << std::endl;
        exit(-1);
    }

    input_file.open(argv[1], std::ios::in);

    if (!input_file.good()) {
        std::cout << "error: could not read <input_file>" << std::endl;
        exit(-1);
    }

    std::string output_file_name = argv[2];
    if (std::filesystem::exists(output_file_name)) std::filesystem::remove(output_file_name);
    output_file.open(output_file_name, std::ios::in | std::ios::out | std::ios::app);

    if (!output_file.good()) {
        std::cout << "error: could not write to <output_file>" << std::endl;
        exit(-1);
    }

    input_file.read((char*) &n, 5);
    std::cout << "decoding (" << format_size(n) << ")" << std::flush;
    auto t1 = now();

    if (n <= std::numeric_limits<uint32_t>::max()) {
        decode<uint32_t>();
    } else {
        decode<uint64_t>();
    }

    input_file.close();
    output_file.close();
    auto t2 = now();
    log_runtime(t1, t2);
    std::cout << "throughput: " << format_throughput(n, time_diff_ns(t1, t2)) << std::endl;
    std::cout << "peak memory consumption: " << format_size(malloc_count_peak()) << std::endl;
    return 0;
}