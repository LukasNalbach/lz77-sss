#include <fstream>
#include <lz77_sss/lz77_sss.hpp>

int arg_idx = 1;
bool decompress = false;
uint64_t bytes_input;
std::string input;
std::string input_file_path;
std::string output_file_path;
std::string tmp_file_path;
std::ifstream input_file;
std::string encoder = "xz";
std::string decoder;
uint32_t encoding_quality = 2;
uint64_t bytes_compressed;
uint16_t num_threads;

void help(std::string message) {
    if (message != "") std::cout << message << std::endl;
    std::cout << "usage: ssszip [...] <input_file>" << std::endl;
    std::cout << " -d                decompress <input_file> (with extension .ssszip.<encoder>)" << std::endl;
    std::cout << " -o <output_file>  output file path (default: <input_file>.ssszip.<encoder>)" << std::endl;
    std::cout << " -t <threads>      number of threads to use (default: all)" << std::endl;
    std::cout << " -e <encoder>      name of the encoder binary (default: xz)" << std::endl;
    std::cout << " -1 .. -9          encoding quality (default: 2)" << std::endl;
    std::cout << " -h                show help" << std::endl;
    exit(-1);
}

void parse_args(char** argv, int argc) {
    std::string arg = argv[arg_idx++];

    if (arg == "-d") {
        decompress = true;
    } else if (arg == "-o") {
        if (arg_idx >= argc - 1)
            help("error: missing parameter after -o option");
        output_file_path = argv[arg_idx++];
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
    } else if (arg.length() == 2 && arg[0] == '-') {
        encoding_quality = arg[1] - '0';
        if (encoding_quality > 9)
            help("error: invalid encoding quality");
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

    std::function<void(factor)> out_it = [&](factor f){
        tmp_ofile << f;

        if (f.len == 0) {
            tmp_ofile.write(&input[i], f.src);
            i += f.src;
        } else {
            i += f.len;
        }
    };

    lz77_sss<pos_t>::template factorize_approximate<
        skip_phrases, lpf_all_external>(input, out_it,
        {.num_threads = num_threads, .log = true});
}

void encode() {
    auto t1 = now();
    input_file.open(input_file_path);
    if (!input_file.good()) {
        std::cout << "error: could not read <input_file>";
        exit(-1);}
    if (output_file_path == "")
        output_file_path = input_file_path;
    output_file_path += ".ssszip." + encoder;
    input_file.seekg(0, std::ios::end);
    bytes_input = input_file.tellg();
    input_file.seekg(0, std::ios::beg);
    auto t0 = now();
    std::cout << "reading input (" << format_size(bytes_input)
        << ")" << std::flush;
    no_init_resize(input, bytes_input);
    read_from_file(input_file, input.data(), bytes_input);
    input_file.close();
    log_runtime(t0);
    if (bytes_input <= std::numeric_limits<uint32_t>::max())
         encode_gapped<uint32_t>();
    else encode_gapped<uint64_t>();
    if (std::filesystem::exists(output_file_path))
        std::filesystem::remove(output_file_path);
    std::cout << "compressing gapped factorization" << std::flush;
    auto t2 = now();
    std::string cmd = encoder + " -q -c" +
        " -" + std::to_string(encoding_quality) +
        (encoder == "xz" ? " -T " + std::to_string(num_threads) : "") +
        " " + tmp_file_path + " > " + output_file_path;
    system(cmd.c_str());
    auto t3 = now();
    log_runtime(t2, t3);
    uint64_t time_total = time_diff_ns(t1, t3);
    uint64_t bytes_compressed = std::filesystem::file_size(output_file_path);
    std::cout << "total time: " << format_time(time_total) << std::endl;
    std::cout << "throughput: " << format_throughput(bytes_input, time_total) << std::endl;
    std::cout << "peak memory consumption: " << format_size(malloc_count_peak()) << std::endl;
    std::cout << "output file size: " << format_size(bytes_compressed) << std::endl;
    std::cout << "compression ratio: " << bytes_input / (double) bytes_compressed << std::endl;
}

template <typename pos_t>
void decode_gapped(std::fstream& tmp_input_file) {
    auto t1 = now();
    using factor = lz77_sss<pos_t>::factor;
    std::fstream output_file(output_file_path,
        std::ios::in | std::ios::out | std::ios::app);
    std::string buff;
    factor f;
    pos_t pos_input = 9;
    pos_t pos_output = 0;
    tmp_input_file.read((char*)&bytes_input, 8);
    std::cout << "decompressing input (" <<
        format_size(bytes_compressed) << ")" << std::flush;
    
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
    std::cout << std::endl;
    std::cout << "total time: " << format_time(time_total) << std::endl;
    std::cout << "throughput: " << format_throughput(bytes_input, time_total) << std::endl;
    std::cout << "peak memory consumption: " << format_size(malloc_count_peak()) << std::endl;
    std::cout << "output file size: " << format_size(bytes_input) << std::endl;
    std::cout << "compression ratio: " << bytes_input / (double) bytes_compressed << std::endl;
}

void decode() {
    bytes_compressed = std::filesystem::file_size(input_file_path);
    encoder = input_file_path.substr(input_file_path.find_last_of(".") + 1);
    if (input_file_path.substr(input_file_path.length() - encoder.length() - 8, 7) != ".ssszip")
        help("error: <input_file> does not have extension .ssszip.<encoder>");
    system((encoder + " -d -c -q " + input_file_path + " > " + tmp_file_path).c_str());
    if (output_file_path == "")
        output_file_path = input_file_path.substr(
        0, input_file_path.length() - encoder.length() - 8);
    if (std::filesystem::exists(output_file_path))
        std::filesystem::remove(output_file_path);
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
    tmp_file_path = std::filesystem::temp_directory_path().string()
        + "/tmp_" + random_alphanumeric_string(10);
    if (decompress) decode(); else encode();
}