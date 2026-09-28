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
#include <lz77_sss/inst.hpp>

std::string alz_binary()
{
    std::error_code ec;
    std::filesystem::path exe = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (ec) return "alz";
    return (exe.parent_path() / "alz").string();
}

uint64_t bytes_input;
std::string input_file_path;
std::string log_file_path;
std::string tmp_file_path;
std::ifstream input_file;
uint16_t min_threads;
uint16_t max_threads;

void help(std::string message)
{
    if (message != "") std::cout << message << std::endl;
    std::cout << "usage: zip-bench [-m <m_file>] <input_file> <min_threads> <max_threads>" << std::endl;
    std::cout << " -m <m_file>    append results to <m_file>" << std::endl;
    std::cout << " <min_threads>  smallest thread count for the multithreaded compressors" << std::endl;
    std::cout << " <max_threads>  largest thread count; the count doubles from <min_threads>" << std::endl;
    exit(-1);
}

uint64_t child_peak_rss_kib()
{
    std::ifstream log_file(log_file_path);
    std::string log_file_str;
    uint64_t log_file_length = std::filesystem::file_size(log_file_path);
    no_init_resize(log_file_str, log_file_length);
    log_file.read(log_file_str.data(), log_file_length);
    log_file.close();
    std::filesystem::remove(log_file_path);
    std::string str_to_find = "Maximum resident set size (kbytes): ";
    uint64_t pos = log_file_str.find(str_to_find);
    if (pos == std::string::npos) return 0;
    uint64_t beg = pos + str_to_find.length();
    uint64_t len = log_file_str.find("\n", beg) - beg;
    return atol(log_file_str.substr(beg, len).c_str());
}

void bench(std::string compressor, bool use_multiple_threads, uint32_t param = 0)
{
    std::string compressor_label = compressor + (param != 0 ? ("_" + std::to_string(param)) : "");
    std::string binary = compressor == "alz" ? alz_binary() : compressor;

    if (system(("command -v " + binary + " > /dev/null 2>&1").c_str()) != 0) {
        std::cout << "skipping " << compressor_label << ": " << binary << " not found" << std::endl;
        return;
    }

    uint32_t min_threads_local = use_multiple_threads ? min_threads : 1;
    uint32_t max_threads_local = use_multiple_threads ? max_threads : 1;

    for (uint16_t num_threads = min_threads_local; num_threads <= max_threads_local; num_threads *= 2) {
        std::cout << "benchmarking " << compressor_label << " using " << num_threads << " threads" << std::flush;
        std::string output_file_path = input_file_path + ".compressed";
        std::filesystem::remove(output_file_path);
        std::string num_thr_str = std::to_string(num_threads);
        std::string cpu_list;

        for (uint32_t i = 0; i < num_threads; i++)
            if (num_threads == omp_get_max_threads())
                cpu_list += std::to_string(i) + ",";
            else
                cpu_list += std::to_string(2 * i) + ",";
        cpu_list.resize(cpu_list.length() - 1);
        std::string cmd = "/usr/bin/time -v taskset -c " + cpu_list + " " +
            (compressor == "alz" ? "env OMP_NUM_THREADS=" + num_thr_str + " " : "") + binary +
            (compressor == "alz" ? (
                " " + input_file_path + " -o " + output_file_path +
                " -s " + std::to_string(param) + " > /dev/null"
            ) : (compressor == "7z" ? (
                " a -m0=lzma2 -mmt" + num_thr_str + " " +
                output_file_path + " " + input_file_path + " > /dev/null"
            ) : (compressor == "bsc" ? (
                " e " + input_file_path + " " + output_file_path + " -b" +
                std::to_string(param) + " -e2 > /dev/null"
            ) : (
                " -k -c -q" +
                (compressor == "xz" ? " -T " + num_thr_str : "") +
                (compressor == "zstd" ? " -T" + num_thr_str : "") +
                " " + input_file_path + " > " + output_file_path))))
            + " 2> " + log_file_path;
        auto t1 = now();
        if (system(cmd.c_str())) {}
        auto t2 = now();

        uint64_t memory_peak_compress = child_peak_rss_kib() * 1024;
        uint64_t time_compress = time_diff_ns(t1, t2);

        if (compressor == "alz") {
            std::filesystem::rename(output_file_path, output_file_path + ".alz");
            cmd = "(/usr/bin/time -v taskset -c " + cpu_list + " bsc e " + output_file_path +
                ".alz " + output_file_path + " -b2047 -e2 > /dev/null) 2> " + log_file_path;

            t1 = now();
            if (system(cmd.c_str())) {}
            t2 = now();

            memory_peak_compress = std::max(memory_peak_compress, child_peak_rss_kib() * 1024);
            time_compress += time_diff_ns(t1, t2);
        }

        uint64_t bytes_compressed = std::filesystem::file_size(output_file_path);
        double compression_ratio = bytes_input / (double)bytes_compressed;
        std::cout << std::endl;
        std::cout << "time = " << format_time(time_compress) << std::endl;
        std::cout << "throughput = " << format_throughput(bytes_input, time_compress) << std::endl;
        std::cout << "peak memory consumption = " << format_size(memory_peak_compress) << std::endl;
        std::cout << "output file size = " << format_size(bytes_compressed) << std::endl;
        std::cout << "compression ratio = " << compression_ratio << std::endl;
        std::cout << std::endl;

        if (result_log::out.is_open()) {
            result_log::out << "RESULT"
                << " text_name=" << result_log::text_name
                << " type=compress"
                << " num_threads=" << num_threads
                << " n=" << bytes_input
                << " compressor=" << compressor_label
                << " time=" << time_compress
                << " throughput=" << throughput_mb_per_s(bytes_input, time_compress)
                << " mem_peak=" << memory_peak_compress
                << " bytes_comp=" << bytes_compressed
                << " comp_ratio=" << compression_ratio
                << std::endl;
        }

        if (num_threads != min_threads_local) continue;
        uint64_t memory_peak_decompress = 0;
        uint64_t time_decompress = 0;
        std::cout << "decompressing" << std::flush;

        if (compressor == "alz") {
            std::string cmd2 = "(/usr/bin/time -v taskset -c " + cpu_list + " bsc d " +
                output_file_path + " " + output_file_path + ".alz > /dev/null) 2> " + log_file_path;

            t1 = now();
            if (system(cmd2.c_str())) {}
            t2 = now();

            memory_peak_decompress = child_peak_rss_kib() * 1024;
            time_decompress = time_diff_ns(t1, t2);
        }

        std::string cmd3 = "(/usr/bin/time -v " + binary +
            (compressor == "alz" ?
                (" " + output_file_path + " -d -o " + tmp_file_path + " > /dev/null") :
            (compressor == "7z" ?
                (" e -o/tmp/ " + output_file_path + " " + result_log::text_name + " > /dev/null") :
            (compressor == "bsc" ? (
                (" d " + output_file_path + " " + tmp_file_path + " > /dev/null")
            ) : (" -c -d -q " + output_file_path + " > " + tmp_file_path))))
            + ") 2> " + log_file_path;
        t1 = now();
        if (system(cmd3.c_str())) {}
        t2 = now();

        std::filesystem::remove(tmp_file_path);
        std::filesystem::remove(output_file_path);
        memory_peak_decompress = std::max(memory_peak_decompress, child_peak_rss_kib() * 1024);
        time_decompress += time_diff_ns(t1, t2);
        std::cout << std::endl;
        std::cout << "time = " << format_time(time_decompress) << std::endl;
        std::cout << "throughput = " << format_throughput(bytes_input, time_decompress) << std::endl;
        std::cout << "peak memory consumption = " << format_size(memory_peak_decompress) << std::endl;
        std::cout << std::endl;

        if (result_log::out.is_open()) {
            result_log::out << "RESULT"
                << " text_name=" << result_log::text_name
                << " type=decompress"
                << " n=" << bytes_input
                << " compressor=" << compressor_label
                << " time=" << time_decompress
                << " throughput=" << throughput_mb_per_s(bytes_input, time_decompress)
                << " mem_peak=" << memory_peak_decompress
                << " bytes_comp=" << bytes_compressed
                << " comp_ratio=" << compression_ratio
                << std::endl;
        }
    }
}

int main(int argc, char** argv)
{
    int arg_idx = 1;

    if (arg_idx + 1 < argc && std::string(argv[arg_idx]) == "-m") {
        result_log::path = argv[arg_idx + 1];
        arg_idx += 2;
    }

    if (argc - arg_idx != 3 || argv[arg_idx][0] == '-') help("");
    input_file_path = argv[arg_idx];
    input_file.open(input_file_path);
    if (!input_file.good()) help("error: could not read <input_file>");
    min_threads = atoi(argv[arg_idx + 1]);
    max_threads = atoi(argv[arg_idx + 2]);

    if (min_threads == 0 || max_threads == 0 ||
        min_threads > max_threads ||
        max_threads > omp_get_max_threads()
    ) help("error: invalid number of threads");

    if (result_log::path != "") {
        result_log::out.open(result_log::path, std::ofstream::app);
        if (!result_log::out.good()) help("error: could not write to <m_file>");
    }

    for (std::string tool : { std::string("/usr/bin/time"), std::string("taskset") })
        if (system(("command -v " + tool + " > /dev/null 2>&1").c_str()) != 0)
            help("error: " + tool + " not found (on Ubuntu: apt install time util-linux)");

    bytes_input = std::filesystem::file_size(input_file_path);
    std::cout << "input file size = " << format_size(bytes_input) << std::endl;
    result_log::text_name = input_file_path.substr(input_file_path.find_last_of("/\\") + 1);
    log_file_path = std::filesystem::temp_directory_path().string()
        + "/log_" + random_alphanumeric_string(10);
    tmp_file_path = std::filesystem::temp_directory_path().string() + "/" + result_log::text_name;

    bench("alz", true, 6);
    bench("bsc", true, 2047);
    bench("lz4", false);
    bench("7z", true);
    bench("gzip", false);
    bench("bzip2", false);
    bench("xz", true);
    bench("zstd", true);

    return 0;
}
