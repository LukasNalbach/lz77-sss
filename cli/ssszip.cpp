#include <fstream>
#include <lz77_sss/lz77_sss.hpp>

using time_point_t = std::chrono::steady_clock::time_point;

time_point_t t1, t2, t3;
int arg_idx = 1;
bool decompress = false;
uint64_t bytes_input;
std::string text_name;
std::string input;
std::string input_file_path;
std::string output_file_path;
std::string result_file_path;
std::string tmp_file_path;
std::string log_file_path;
std::ifstream input_file;
std::ofstream result_file;
std::string encoder = "zstd";
uint32_t encoding_quality = 4;
uint64_t bytes_compressed;
uint16_t num_threads;
uint8_t logs = 1;
uint64_t gaps_length = 0;
bool keep_input_file = false;

void help(std::string message) {
    if (logs >= 1) {
        if (message != "") std::cout << message << std::endl;
        std::cout << "usage: ssszip [...] <input_file>" << std::endl;
        std::cout << " -d                decompress <input_file> (with extension .ssszip.<encoder>)" << std::endl;
        std::cout << " -o <output_file>  output file path (default: <input_file>.ssszip.<encoder>)" << std::endl;
        std::cout << " -t <threads>      number of threads to use (default: all)" << std::endl;
        std::cout << " -e <encoder>      name of the encoder binary (default: zstd)" << std::endl;
        std::cout << " -0/-1/-2/...      encoding quality (default: 4)" << std::endl;
        std::cout << " -k                keep (don't delete) <input file>" << std::endl;
        std::cout << " -q                quiet mode (disables all logs)" << std::endl;
        std::cout << " -v                shows verbose information" << std::endl;
        std::cout << " -r <result_file>  write results to <result_file>" << std::endl;
        std::cout << " -h                show help" << std::endl;
    }
    exit(-1);
}

void parse_args(char** argv, int argc) {
    std::string arg = argv[arg_idx++];

    if (arg == "-k") {
        keep_input_file = true;
    } else if (arg == "-q") {
        if (logs != 1)
            help("error: -q and -v cannot be combined");
        logs = 0;
    } else if (arg == "-v") {
        if (logs != 1)
            help("error: -q and -v cannot be combined");
        logs = 2;
    } else if (arg == "-d") {
        decompress = true;
    } else if (arg == "-o") {
        if (arg_idx >= argc - 1)
            help("error: missing parameter after -o option");
        output_file_path = argv[arg_idx++];
    } else if (arg == "-r") {
        if (arg_idx >= argc - 1)
            help("error: missing parameter after -r option");
        result_file_path = argv[arg_idx++];
    } else if (arg == "-t") {
        if (arg_idx >= argc - 1)
            help("error: missing parameter after -t option");
        num_threads = std::max<uint16_t>(1, atoi(argv[arg_idx++]));
        if (num_threads > omp_get_max_threads())
            help("error: requested too many threads");
    } else if (arg == "-e") {
        if (arg_idx >= argc - 1)
            help("error: missing parameter after -e option");
        encoder = argv[arg_idx++];
    } else if (arg == "-h") {
        help("");
    } else if (2 <= arg.length() && arg.length() <= 3 && arg[0] == '-' &&
        ('0' <= arg[1] && arg[1] <= '9') && (arg.length() == 2 || ('0' <= arg[2] && arg[2] <= '9'))
    ) {
        encoding_quality = arg[1] - '0';
        if (arg.length() == 3)
            encoding_quality = 10 * encoding_quality + arg[2] - '0';
    } else {
        help("error: unrecognized '" + arg + "' option");
    }
}

template <typename pos_t>
void encode_gapped() {
    using factor = lz77_sss<pos_t>::factor;
    std::ofstream tmp_ofile(tmp_file_path);
    bool is_64_bit = std::is_same_v<pos_t, uint64_t>;
    tmp_ofile.write((char*)&is_64_bit, 1);
    tmp_ofile.write((char*)&bytes_input, 8);
    pos_t i = 0;
    bool gap = false;
    pos_t gap_beg = 0;
    factor gap_lst {.src = 0, .len = 0};

    std::function<void(factor)> out_it = [&](factor f){
        if (f.len == 0) {
            if (gap) {
                gap_lst.src += f.src;
            } else {
                gap_lst = f;
                gap_beg = i;
            }

            i += f.src;
            gap = true;
        } else if (gap && f.len <= sizeof(factor)) {
            gap_lst.src += f.len;
            i += f.len;
            gap = true;
        } else {
            if (gap) {
                tmp_ofile << gap_lst;
                tmp_ofile.write(&input[gap_beg], gap_lst.src);
                gaps_length += gap_lst.src;
            }

            tmp_ofile << f;
            i += f.len;
            gap = false;
        }
    };

    lz77_sss<pos_t>::template factorize_approximate<
        skip_phrases, lpf_all_external>(input, out_it,
        {.num_threads = num_threads, .log = logs == 2});
    
    if (gap) {
        tmp_ofile << gap_lst;
        tmp_ofile.write(&input[gap_beg], gap_lst.src);
        gaps_length += gap_lst.src;
    }
}

uint64_t peak_memory_usage() {
    std::ifstream log_file(log_file_path);
    std::string log_file_str;
    uint64_t log_file_length = std::filesystem::file_size(log_file_path);
    no_init_resize(log_file_str, log_file_length);
    log_file.read(log_file_str.data(), log_file_length);
    log_file.close();
    std::string str_to_find = "Maximum resident set size (kbytes): ";
    uint64_t beg = log_file_str.find(str_to_find) + str_to_find.length();
    uint64_t len = log_file_str.find("\n", beg) - beg;
    return atol(log_file_str.substr(beg, len).c_str());
}

void encode() {
    t1 = now();
    input_file.open(input_file_path);
    if (!input_file.good())
        help("error: could not read <input_file>");
    if (output_file_path == "")
        output_file_path = input_file_path;
    output_file_path += ".ssszip." + encoder;
    if (std::filesystem::exists(output_file_path))
        help("error: <output_file> already exists");
    input_file.seekg(0, std::ios::end);
    bytes_input = input_file.tellg();
    input_file.seekg(0, std::ios::beg);
    if (logs == 2) std::cout << "reading input (" <<
        format_size(bytes_input) << ")" << std::flush;
    no_init_resize(input, bytes_input);
    read_from_file(input_file, input.data(), bytes_input);
    input_file.close();
    if (logs == 2) log_runtime(t1);
    if (bytes_input <= std::numeric_limits<uint32_t>::max())
         encode_gapped<uint32_t>();
    else encode_gapped<uint64_t>();
    input.clear();
    input.shrink_to_fit();
    if (logs == 2) {
        uint64_t bytes_gapped = std::filesystem::file_size(tmp_file_path);
        std::cout << "size: " << format_size(bytes_gapped);
        std::cout << ", relative length of the gaps: " << 100.0 *
            (gaps_length / (double) bytes_input) << " %" << std::endl;
        std::cout << "compressing gapped factorization" << std::flush;
    }
    t2 = now();
    std::string cmd = "(/usr/bin/time -v " + encoder + " -c" +
        (logs <= 1 ? " -q" : " -v") +
        " -" + std::to_string(encoding_quality) +
        (encoder == "xz" ? " -T " + std::to_string(num_threads) : "") +
        (encoder == "zstd" ? " -T" + std::to_string(num_threads) : "") +
        " " + tmp_file_path + " > " + output_file_path + ") 2> " + log_file_path;
    system(cmd.c_str());
    t3 = now();
    uint64_t encoding_peak = peak_memory_usage() * 1000;
    uint64_t memory_peak = std::max(encoding_peak, malloc_count_peak());
    uint64_t time_total = time_diff_ns(t1, t3);
    uint64_t bytes_compressed = std::filesystem::file_size(output_file_path);
    double compression_ratio = bytes_input / (double) bytes_compressed;
    std::filesystem::remove(log_file_path);
    
    if (logs == 2) {
        std::cout << ", in " << format_time(time_diff_ns(t2, t3)) << std::endl;
        std::cout << "peak memory usage: " << format_size(encoding_peak) << std::endl;
        std::cout << "total time: " << format_time(time_total) << std::endl;
        std::cout << "total throughput: " << format_throughput(bytes_input, time_total) << std::endl;
        std::cout << "total peak memory consumption: " << format_size(memory_peak) << std::endl;
        std::cout << "output file size: " << format_size(bytes_compressed) << std::endl;
        std::cout << "compression ratio: " << compression_ratio << std::endl;
    }

    if (result_file_path != "") {
        result_file.open(result_file_path, std::ios_base::app);
        text_name = input_file_path.substr(input_file_path.find_last_of("/\\") + 1);

        result_file << "RESULT"
            << " text_name=" << text_name
            << " type=compress"
            << " num_threads=" << num_threads
            << " n=" << bytes_input
            << " encoder=ssszip_" << encoder
            << " time=" << time_total
            << " throughput=" << throughput(bytes_input, time_total)
            << " mem_peak=" << memory_peak
            << " bytes_comp=" << bytes_compressed
            << " comp_ratio=" << compression_ratio
            << std::endl;
    }
}

template <typename pos_t>
void decode_gapped(std::fstream& tmp_input_file) {
    t2 = now();
    using factor = lz77_sss<pos_t>::factor;
    std::fstream output_file(output_file_path,
        std::ios::in | std::ios::out | std::ios::app);
    std::string buff;
    factor f;
    pos_t pos_input = 9;
    pos_t pos_output = 0;
    tmp_input_file.read((char*)&bytes_input, 8);
    uint64_t bytes_gapped = std::filesystem::file_size(tmp_file_path);
    if (logs == 2) std::cout << "reverting gapped factorization (" <<
        format_size(bytes_gapped) << ")" << std::flush;
    
    while (pos_output < bytes_input) {
        tmp_input_file >> f;
        pos_input += sizeof(factor);

        if (f.len == 0) {
            copy_buffered(
                tmp_input_file, output_file,
                buff, pos_input, pos_output, f.src);
            pos_input += f.src;
            pos_output += f.src;
        } else {
            copy_buffered(
                output_file, output_file,
                buff, f.src, pos_output, f.len);
            pos_output += f.len;
        }
    }

    uint64_t time_total = time_diff_ns(t1, now());
    uint64_t decoding_peak = peak_memory_usage() * 1000;
    uint64_t memory_peak = std::max(decoding_peak, malloc_count_peak());
    uint64_t compression_ratio = bytes_input / (double) bytes_compressed;

    if (logs == 2) {
        log_runtime(t2);
        std::cout << "total time: " << format_time(time_total) << std::endl;
        std::cout << "throughput: " << format_throughput(bytes_input, time_total) << std::endl;
        std::cout << "peak memory consumption: " << format_size(memory_peak) << std::endl;
        std::cout << "output file size: " << format_size(bytes_input) << std::endl;
        std::cout << "compression ratio: " << compression_ratio << std::endl;
    }

    if (result_file_path != "") {
        result_file.open(result_file_path, std::ios_base::app);
        text_name = output_file_path.substr(output_file_path.find_last_of("/\\") + 1);

        result_file << "RESULT"
            << " text_name=" << text_name
            << " type=decompress"
            << " n=" << bytes_input
            << " encoder=ssszip_" << encoder
            << " time=" << time_total
            << " throughput=" << throughput(bytes_input, time_total)
            << " mem_peak=" << memory_peak
            << " bytes_comp=" << bytes_compressed
            << " comp_ratio=" << compression_ratio
            << std::endl;
    }
}

void decode() {
    t1 = now();
    if (output_file_path == "")
        output_file_path = input_file_path.substr(
        0, input_file_path.length() - encoder.length() - 8);
    if (std::filesystem::exists(output_file_path))
        help("error: <output_file> already exists");
    bytes_compressed = std::filesystem::file_size(input_file_path);
    encoder = input_file_path.substr(input_file_path.find_last_of(".") + 1);
    if (input_file_path.substr(input_file_path.length() - encoder.length() - 8, 7) != ".ssszip")
        help("error: <input_file> does not have extension .ssszip.<encoder>");
    if (logs == 2) std::cout << "decompressing gapped factorization (" <<
        format_size(bytes_compressed) << ")" << std::flush;
    std::string cmd = "(/usr/bin/time -v " + encoder + " -d -c -q " +
        input_file_path + " > " + tmp_file_path + ") 2> " + log_file_path;
    system(cmd.c_str());
    if (logs == 2) log_runtime(t1);
    std::fstream tmp_input_file(tmp_file_path, std::ios::in);
    bool is_64_bit;
    tmp_input_file.read((char*)&is_64_bit, 1);
    if (is_64_bit) decode_gapped<uint64_t>(tmp_input_file);
              else decode_gapped<uint32_t>(tmp_input_file);
    tmp_input_file.close();
    std::filesystem::remove(tmp_file_path);
}

int main(int argc, char** argv) {
    num_threads = omp_get_max_threads();
    while (arg_idx < argc - 1) parse_args(argv, argc);
    input_file_path = argv[arg_idx];
    if (!std::filesystem::exists(input_file_path)) help("error: <input_file> does not exist");
    tmp_file_path = std::filesystem::temp_directory_path().string()
        + "/tmp_" + random_alphanumeric_string(10);
    log_file_path = std::filesystem::temp_directory_path().string()
        + "/log_" + random_alphanumeric_string(10);
    if (decompress) decode(); else encode();
    if (!keep_input_file) std::filesystem::remove(input_file_path);
}