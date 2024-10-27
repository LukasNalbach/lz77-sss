#include <fstream>
#include <lz77_sss/lz77_sss.hpp>

uint64_t bytes_input;
std::string input_file_path;
std::string log_file_path;
std::string tmp_file_path;
std::string result_file_path;
std::ifstream input_file;
std::string text_name;
uint16_t min_threads;
uint16_t max_threads;

void help(std::string message)
{
    if (message != "") std::cout << message << std::endl;
    std::cout << "usage: zip-bench <input_file> <min_threads> <max_threads>" << std::endl;
    exit(-1);
}

uint64_t peak_memory_usage()
{
    std::ifstream log_file(log_file_path);
    std::string log_file_str;
    uint64_t log_file_length = std::filesystem::file_size(log_file_path);
    no_init_resize(log_file_str, log_file_length);
    log_file.read(log_file_str.data(), log_file_length);
    log_file.close();
    std::filesystem::remove(log_file_path);
    std::string str_to_find = "Maximum resident set size (kbytes): ";
    uint64_t beg = log_file_str.find(str_to_find) + str_to_find.length();
    uint64_t len = log_file_str.find("\n", beg) - beg;
    return atol(log_file_str.substr(beg, len).c_str());
}

void bench(std::string encoder, bool use_multiple_threads)
{
    uint32_t min_threads_local = use_multiple_threads ? min_threads : 1;
    uint32_t max_threads_local = use_multiple_threads ? max_threads : 1;

    for (uint16_t num_threads = min_threads_local; num_threads <= max_threads_local; num_threads *= 2) {
        std::cout << "benchmarking " << encoder << " using " << num_threads << " threads" << std::flush;
        std::string output_file_path = input_file_path + ".compressed";
        std::filesystem::remove(output_file_path);
        std::string cpu_list;
        for (uint32_t i = 0; i < num_threads; i++)
            if (num_threads == omp_get_max_threads())
                cpu_list += std::to_string(i) + ",";
            else
                cpu_list += std::to_string(2 * i) + ",";
        cpu_list.resize(cpu_list.length() - 1);
        std::string cmd = "(/usr/bin/time -v taskset -c " + cpu_list + " " + encoder +
            (encoder == "7z" ? (
                " a -m0=lzma2 -mmt" + std::to_string(num_threads) + " " +
                output_file_path + " " + input_file_path + " > /dev/null"
            ) : (
                " -k -c -q" +
                (encoder == "xz" ? " -T " + std::to_string(num_threads) : "") +
                (encoder == "zstd" ? " -T" + std::to_string(num_threads) : "") +
                " " + input_file_path + " > " + output_file_path))
            + ") 2> " + log_file_path;
        auto t1 = now();
        system(cmd.c_str());
        auto t2 = now();

        uint64_t memory_peak_compress = peak_memory_usage() * 1000;
        uint64_t time_compress = time_diff_ns(t1, t2);
        uint64_t bytes_compressed = std::filesystem::file_size(output_file_path);
        double compression_ratio = bytes_input / (double)bytes_compressed;
        std::cout << std::endl;
        std::cout << "time: " << format_time(time_compress) << std::endl;
        std::cout << "throughput: " << format_throughput(bytes_input, time_compress) << std::endl;
        std::cout << "peak memory consumption: " << format_size(memory_peak_compress) << std::endl;
        std::cout << "output file size: " << format_size(bytes_compressed) << std::endl;
        std::cout << "compression ratio: " << compression_ratio << std::endl;
        std::cout << std::endl;

        if (result_file_path != "") {
            std::ofstream result_file(result_file_path, std::ofstream::app);

            result_file << "RESULT"
                << " text_name=" << text_name
                << " type=compress"
                << " num_threads=" << num_threads
                << " n=" << bytes_input
                << " encoder=" << encoder
                << " time=" << time_compress
                << " throughput=" << throughput(bytes_input, time_compress)
                << " mem_peak=" << memory_peak_compress
                << " bytes_comp=" << bytes_compressed
                << " comp_ratio=" << compression_ratio
                << std::endl;
        }

        if (num_threads != 1) continue;
        std::cout << "decompressing" << std::flush;
        std::string cmd2 = "(/usr/bin/time -v " + encoder +
            (encoder == "7z" ?
                (" e -o/tmp/ " + output_file_path + " " + text_name + " > /dev/null") :
                (" -c -d -q " + output_file_path + " > " + tmp_file_path))
            + ") 2> " + log_file_path;
        t1 = now();
        system(cmd2.c_str());
        t2 = now();

        std::filesystem::remove(tmp_file_path);
        std::filesystem::remove(output_file_path);
        uint64_t memory_peak_decompress = peak_memory_usage() * 1000;
        uint64_t time_decompress = time_diff_ns(t1, t2);
        std::cout << "time: " << format_time(time_decompress) << std::endl;
        std::cout << "throughput: " << format_throughput(bytes_input, time_decompress) << std::endl;
        std::cout << "peak memory consumption: " << format_size(memory_peak_decompress) << std::endl;
        std::cout << std::endl;

        if (result_file_path != "") {
            std::ofstream result_file(result_file_path, std::ofstream::app);
            
            result_file << "RESULT"
                << " text_name=" << text_name
                << " type=decompress"
                << " n=" << bytes_input
                << " encoder=" << encoder
                << " time=" << time_decompress
                << " throughput=" << throughput(bytes_input, time_decompress)
                << " mem_peak=" << memory_peak_decompress
                << " bytes_comp=" << bytes_compressed
                << " comp_ratio=" << compression_ratio
                << std::endl;
        }
    }
}

int main(int argc, char** argv)
{
    if (!(4 <= argc && argc <= 5))
        help("error: invalid number of parameters");

    input_file_path = argv[1];
    input_file.open(input_file_path);

    if (!input_file.good())
        help("error: could not read <file>");

    min_threads = atoi(argv[2]);
    max_threads = atoi(argv[3]);

    if (min_threads == 0 || max_threads == 0 ||
        min_threads > max_threads ||
        max_threads > omp_get_max_threads()
    ) help("error: invalid number of threads");

    if (argc == 5)
        result_file_path = argv[4];

    bytes_input = std::filesystem::file_size(input_file_path);
    std::cout << "input_file size: " << format_size(bytes_input) << std::endl;
    text_name = input_file_path.substr(input_file_path.find_last_of("/\\") + 1);
    log_file_path = std::filesystem::temp_directory_path().string()
        + "/log_" + random_alphanumeric_string(10);
    tmp_file_path = std::filesystem::temp_directory_path().string() + "/" + text_name;

    bench("lz4", false);
    bench("7z", true);
    bench("gzip", false);
    bench("bzip2", false);
    bench("xz", true);
    bench("zstd", true);
}