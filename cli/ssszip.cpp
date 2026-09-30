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

#include <algorithm>
#include <array>
#include <bit>
#include <csignal>
#include <cstdio>
#include <fstream>
#include <numeric>

#include <lz77_sss/inst.hpp>
#include <lz77_sss/misc/vbyte.hpp>
#include <lz77_sss/misc/text_reader.hpp>
#include <lz77_sss/misc/file_decoder.hpp>

using time_point_t = std::chrono::steady_clock::time_point;
static constexpr uint64_t min_lpf_len = 128;
static constexpr double min_two_bit_share = 0.9;
static constexpr uint8_t gapped_fasta_flag = 2;
static constexpr uint8_t gapped_two_bit_flag = 4;
static constexpr uint8_t gapped_streams_flag = 8;
static constexpr uint64_t lit_write_buffer = uint64_t { 1 } << 22;
static constexpr uint64_t token_write_buffer = uint64_t { 1 } << 16;
static constexpr uint64_t stream_read_buffer = uint64_t { 1 } << 20;

time_point_t time_start;
int arg_idx = 1;
bool decompress_mode = false;
uint64_t bytes_input;
std::string input_file_path;
std::string output_file_path;
std::string tmp_file_path;
std::string log_file_path;
std::ifstream input_file;
std::string postcompressor = "bsc";
uint32_t post_compression_quality = 4;
bool post_compression_quality_given = false;
static constexpr double bsc_bytes_per_block_byte = 5.5;
static constexpr uint64_t bsc_min_block = 10000;
static constexpr uint64_t bsc_max_block = uint64_t { 2047 } << 20;
static constexpr uint64_t feed_block = uint64_t { 1 } << 20;
uint64_t bytes_compressed;
uint16_t num_threads;
bool quiet = false;
bool verbose_given = false;
uint64_t gaps_length = 0;
bool two_bit_literals = false;
bool decode_in_ram = false;
text_encoding encoding = auto_encoding;
fasta_mode fasta = fasta_auto;
std::string err_file_path;

enum class postcompressor_kind {
    zstd, xz, lzma, gzip, pigz, bzip2, pbzip2, lbzip2, brotli, lz4, lzop, lzip, plzip, bzip3, sevenzip, bsc
};

struct postcompressor_spec {
    postcompressor_kind kind;
    const char* name;
    const char* binary;
    uint32_t min_quality;
    uint32_t max_quality;
    uint32_t default_quality;
};

static constexpr std::array<postcompressor_spec, 16> postcompressors { {
    { postcompressor_kind::bsc, "bsc", "bsc", 1, 2047, 0 },
    { postcompressor_kind::zstd, "zstd", "zstd", 1, 22, 4 },
    { postcompressor_kind::xz, "xz", "xz", 0, 9, 6 },
    { postcompressor_kind::lzma, "lzma", "xz", 0, 9, 6 },
    { postcompressor_kind::gzip, "gzip", "gzip", 1, 9, 6 },
    { postcompressor_kind::pigz, "pigz", "pigz", 1, 9, 6 },
    { postcompressor_kind::bzip2, "bzip2", "bzip2", 1, 9, 9 },
    { postcompressor_kind::pbzip2, "pbzip2", "pbzip2", 1, 9, 9 },
    { postcompressor_kind::lbzip2, "lbzip2", "lbzip2", 1, 9, 9 },
    { postcompressor_kind::brotli, "brotli", "brotli", 0, 11, 11 },
    { postcompressor_kind::lz4, "lz4", "lz4", 1, 12, 1 },
    { postcompressor_kind::lzop, "lzop", "lzop", 1, 9, 3 },
    { postcompressor_kind::lzip, "lzip", "lzip", 0, 9, 6 },
    { postcompressor_kind::plzip, "plzip", "plzip", 0, 9, 6 },
    { postcompressor_kind::bzip3, "bzip3", "bzip3", 1, 511, 16 },
    { postcompressor_kind::sevenzip, "7z", "7z", 0, 9, 5 }
} };

const postcompressor_spec* spec = nullptr;

enum class progress_source {
    pipe_input,
    tool_stdout,
    tool_stderr
};

struct postcompressor_command {
    std::string command;
    progress_source progress;
};

const postcompressor_spec* find_postcompressor(const std::string& name)
{
    for (const postcompressor_spec& pc : postcompressors) {
        if (name == pc.name) return &pc;
    }

    return nullptr;
}

bool binary_installed(const std::string& name)
{
    const char* path = std::getenv("PATH");
    if (path == nullptr) return false;
    #ifdef _WIN32
    const char separator = ';';
    const std::string suffix = ".exe";
    #else
    const char separator = ':';
    const std::string suffix;
    #endif
    const std::string dirs(path);

    for (size_t beg = 0; beg <= dirs.size();) {
        size_t end = dirs.find(separator, beg);
        if (end == std::string::npos) end = dirs.size();
        const std::filesystem::path dir = dirs.substr(beg, end - beg);
        std::error_code ec;
        if (!dir.empty() && std::filesystem::is_regular_file(dir / name, ec)) return true;
        if (!dir.empty() && !suffix.empty() && std::filesystem::is_regular_file(dir / (name + suffix), ec)) return true;
        beg = end + 1;
    }

    return false;
}

std::string quote(const std::string& path)
{
    #ifdef _WIN32
    return "\"" + path + "\"";
    #else
    std::string quoted = "'";

    for (char c : path) {
        if (c == '\'') quoted += "'\\''";
        else quoted += c;
    }

    return quoted + "'";
    #endif
}

FILE* open_pipe(const std::string& command, bool write)
{
    #ifdef _WIN32
    return _popen(command.c_str(), write ? "wb" : "rb");
    #else
    return popen(command.c_str(), write ? "w" : "r");
    #endif
}

int close_pipe(FILE* pipe)
{
    #ifdef _WIN32
    return _pclose(pipe);
    #else
    return pclose(pipe);
    #endif
}

void help(std::string message)
{
    if (!quiet) {
        if (message != "") std::cout << message << std::endl;
        std::cout << "usage: ssszip [options] <input_file>" << std::endl;
        std::cout << " -d                    decompress <input_file> (*.ssszip.<postcompressor>)" << std::endl;
        std::cout << " -ram                  keep the whole output in memory while decompressing" << std::endl;
        std::cout << "                       (default: only the last 64 MiB)" << std::endl;
        std::cout << " -o <base_name>        write <base_name>.ssszip.<postcompressor> (default:" << std::endl;
        std::cout << "                       <input_file>); with -d: the output file (default:" << std::endl;
        std::cout << "                       <input_file> without .ssszip.<postcompressor>)" << std::endl;
        std::cout << " -t <threads>          number of threads (default: all)" << std::endl;
        std::string line = " -pc <postcompressor>  postcompressor:";

        for (uint64_t k = 0; k < postcompressors.size(); k++) {
            const std::string name = std::string(postcompressors[k].name) + (k == 0 ? " (default)" : "")
                + (k + 1 < postcompressors.size() ? "," : "");

            if (line.size() + name.size() + 1 > 80) {
                std::cout << line << std::endl;
                line = std::string(23, ' ');
            } else if (line.back() != ' ') {
                line += ' ';
            }

            line += name;
        }

        std::cout << line << std::endl;
        std::cout << " -<quality>            post-compression quality, e.g. -9: higher gives smaller" << std::endl;
        std::cout << "                       output but takes longer (range and default depend on" << std::endl;
        std::cout << "                       the postcompressor); for bsc and bzip3 the block size" << std::endl;
        std::cout << "                       in MB (bsc default: largest block that needs no more" << std::endl;
        std::cout << "                       memory than the factorization, at most 2047)" << std::endl;
        std::cout << " -enc <encoding>       how the text is kept in memory (default: auto):" << std::endl;
        std::cout << "                       plain   one byte per character (fastest)" << std::endl;
        std::cout << "                       packed  fewer bits per character for small alphabets" << std::endl;
        std::cout << "                       split   fewer bits for frequent characters (slowest)" << std::endl;
        std::cout << "                       auto    packed or split if that saves memory, else plain" << std::endl;
        std::cout << " -fasta <mode>         store the header lines of FASTA files separately from" << std::endl;
        std::cout << "                       the sequences: on, off or auto (default: auto, on for" << std::endl;
        std::cout << "                       FASTA files)" << std::endl;
        std::cout << " -q                    print nothing" << std::endl;
        std::cout << " -v                    print progress and statistics (default)" << std::endl;
        std::cout << " -m <m_file>           append results to <m_file>" << std::endl;
        std::cout << " -h                    show help" << std::endl;
    }

    exit(-1);
}

void parse_arg(int argc, char** argv)
{
    std::string arg = argv[arg_idx++];

    if (arg == "-fasta") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -fasta option");
        fasta = parse_fasta_mode(argv[arg_idx++], help);
    } else if (arg == "-q") {
        if (verbose_given) help("error: -q and -v cannot be combined");
        quiet = true;
    } else if (arg == "-v") {
        if (quiet) help("error: -q and -v cannot be combined");
        verbose_given = true;
    } else if (arg == "-d") {
        decompress_mode = true;
    } else if (arg == "-ram") {
        decode_in_ram = true;
    } else if (arg == "-o") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -o option");
        output_file_path = argv[arg_idx++];
    } else if (arg == "-m") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -m option");
        result_log::path = argv[arg_idx++];
    } else if (arg == "-t") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -t option");
        num_threads = std::max<uint16_t>(1, atoi(argv[arg_idx++]));
        if (num_threads > omp_get_max_threads()) help("error: requested too many threads");
    } else if (arg == "-enc") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -enc option");
        encoding = parse_text_encoding(argv[arg_idx++], help);
    } else if (arg == "-pc") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -pc option");
        postcompressor = argv[arg_idx++];
    } else if (arg == "-h") {
        help("");
    } else if (arg[0] == '-') {

        for (uint32_t i = 1; i < arg.size(); i++) {
            if (!std::isdigit(arg[i]))
                help("error: unrecognized '" + arg + "' option");
        }

        post_compression_quality = std::stoi(arg.substr(1));
        post_compression_quality_given = true;
    } else {
        help("error: unrecognized '" + arg + "' option");
    }
}

enum gapped_stream {
    gapped_lits,
    gapped_excs,
    gapped_lit_lens,
    gapped_copy_lens,
    gapped_classes,
    gapped_mantissas,
    num_gapped_streams
};

struct gapped_recent_dists {
    std::array<uint64_t, 3> dist { };

    int find(uint64_t d) const
    {
        for (int k = 0; k < 3; k++) {
            if (dist[k] == d) return k;
        }

        return -1;
    }

    void use(uint64_t d)
    {
        const int k = find(d);
        for (int j = k < 0 ? 2 : k; j > 0; j--) dist[j] = dist[j - 1];
        dist[0] = d;
    }
};

struct gapped_patch {
    uint64_t off;
    uint8_t byte;
};

struct gapped_section {
    uint64_t literals = 0;
    uint64_t regular = 0;
    uint64_t copies = 0;
    uint64_t lead = 0;
    uint64_t trail = 0;
    uint64_t lit_len_bytes = 0;
    uint64_t copy_len_bytes = 0;
    uint64_t mantissa_bits = 0;
    uint64_t exc_bytes = 0;
    uint64_t first_exc = 0;
    uint64_t excs_end = 0;
    bool has_excs = false;
    uint8_t num_new = 0;
    std::array<uint64_t, 3> new_dists { };
    gapped_recent_dists local;
    gapped_recent_dists incoming;
    uint64_t lit_beg = 0;
    uint64_t first_lit_len = 0;
    uint64_t lit_len_off = 0;
    uint64_t copy_len_off = 0;
    uint64_t class_off = 0;
    uint64_t mantissa_off = 0;
    uint64_t exc_off = 0;
    uint64_t prev_excs_end = 0;
    std::vector<gapped_patch> patches;
};

class gapped_byte_writer {
public:
    gapped_byte_writer(positional_writer& writer, uint64_t off, uint16_t t, uint64_t capacity)
        : writer(writer)
        , off(off)
        , capacity(capacity)
        , t(t)
    { }

    uint64_t offset() const { return off + fill; }

    void put(uint8_t byte)
    {
        if (fill == buffer.size()) make_room();
        buffer[fill++] = char(byte);
    }

    void put_vbyte(uint64_t x)
    {
        do {
            uint8_t byte = x & 0x7F;
            x >>= 7;
            if (x) byte |= 0x80;
            put(byte);
        } while (x);
    }

    void put_bytes(const char* data, uint64_t len)
    {
        if (len >= capacity) {
            flush();
            writer.write(data, len, off, t);
            off += len;
            return;
        }

        while (len > 0) {
            if (fill == buffer.size()) make_room();
            const uint64_t take = std::min<uint64_t>(len, buffer.size() - fill);
            std::memcpy(buffer.data() + fill, data, take);
            fill += take;
            data += take;
            len -= take;
        }
    }

    template <typename next_t>
    void put_chars(uint64_t len, next_t next)
    {
        while (len > 0) {
            if (fill == buffer.size()) make_room();
            const uint64_t take = std::min<uint64_t>(len, buffer.size() - fill);
            for (uint64_t k = 0; k < take; k++) buffer[fill + k] = char(next());
            fill += take;
            len -= take;
        }
    }

    void flush()
    {
        if (fill == 0) return;
        writer.write(buffer.data(), fill, off, t);
        off += fill;
        fill = 0;
    }

private:
    void make_room()
    {
        flush();
        if (buffer.size() < capacity) no_init_resize(buffer, capacity);
    }

    positional_writer& writer;
    std::string buffer;
    uint64_t off;
    uint64_t fill = 0;
    uint64_t capacity;
    uint16_t t;
};

class gapped_bit_writer {
public:
    gapped_bit_writer(positional_writer& writer, uint64_t base, uint64_t bit_off, uint16_t t, uint64_t capacity,
        std::vector<gapped_patch>& patches)
        : bytes(writer, base + div_ceil<uint64_t>(bit_off, 8), t, capacity)
        , patches(patches)
        , head_off(base + bit_off / 8)
        , head(bit_off % 8 != 0)
        , bits(bit_off % 8)
    { }

    void put(uint64_t x, uint8_t width)
    {
        if (width > 32) {
            put(x >> 32, width - 32);
            width = 32;
        }

        if (width == 0) return;
        acc = (acc << width) | (x & ((uint64_t { 1 } << width) - 1));
        bits += width;
        used = true;

        while (bits >= 8) {
            bits -= 8;
            const uint8_t byte = uint8_t(acc >> bits);

            if (head) {
                patches.push_back({ head_off, byte });
                head = false;
            } else {
                bytes.put(byte);
            }
        }
    }

    void finish()
    {
        if (used && bits > 0) patches.push_back({ head ? head_off : bytes.offset(), uint8_t(acc << (8 - bits)) });
        bytes.flush();
    }

private:
    gapped_byte_writer bytes;
    std::vector<gapped_patch>& patches;
    uint64_t head_off;
    uint64_t acc = 0;
    bool head;
    bool used = false;
    uint8_t bits;
};

template <typename text_t>
class gapped_encoder {
public:
    gapped_encoder(const text_t& T, const char_histogram& histogram, positional_writer& writer, uint64_t base,
        bool fasta_active)
        : T(T)
        , writer(writer)
        , base(base)
        , fasta_active(fasta_active)
    {
        std::array<uint16_t, 256> order;
        std::iota(order.begin(), order.end(), uint16_t { 0 });
        std::stable_sort(order.begin(), order.end(), [&](uint16_t a, uint16_t b) { return histogram[a] > histogram[b]; });
        uint64_t top = 0;

        for (uint8_t k = 0; k < 4; k++) {
            chars[k] = uint8_t(order[k]);
            top += histogram[order[k]];
        }

        std::sort(chars.begin(), chars.end());

        for (uint8_t k = 0; k < 4; k++) {
            code[chars[k]] = k;
            regular[chars[k]] = true;
        }

        scan = T.size() > 0 && top >= min_two_bit_share * T.size();
    }

    uint64_t num_literals() const { return total_literals; }

    bool uses_two_bit() const { return two_bit; }

    void encode(const lz77_sss::gapped_factorization* view)
    {
        gapped = view;
        sections = std::vector<gapped_section>(view == nullptr ? 0 : view->num_sections());
        const uint64_t num_sections = sections.size();

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
        for (uint64_t s = 0; s < num_sections; s++) encode_section<false>(s, 0);

        arrange();
        const char flags = char(1 | (fasta_active ? gapped_fasta_flag : 0) | (two_bit ? gapped_two_bit_flag : 0) | gapped_streams_flag);
        writer.write(&flags, 1, 0, 0);
        writer.write(header.data(), header.size(), base, 0);

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
        for (uint64_t s = 0; s < num_sections; s++) encode_section<true>(s, omp_get_thread_num());

        std::string last;
        append_vbyte(last, final_lit_len);
        writer.write(last.data(), last.size(), final_lit_len_off, 0);
        std::vector<gapped_patch> patches;

        for (const gapped_section& sect : sections) {
            patches.insert(patches.end(), sect.patches.begin(), sect.patches.end());
        }

        std::sort(patches.begin(), patches.end(), [](const gapped_patch& a, const gapped_patch& b) { return a.off < b.off; });

        for (uint64_t k = 0; k < patches.size();) {
            const uint64_t off = patches[k].off;
            char byte = 0;
            while (k < patches.size() && patches[k].off == off) byte |= char(patches[k++].byte);
            writer.write(&byte, 1, off, 0);
        }
    }

private:
    static uint64_t vbyte_size(uint64_t x)
    {
        return std::max<uint64_t>(1, (std::bit_width(x) + 6) / 7);
    }

    static void append_vbyte(std::string& out, uint64_t x)
    {
        do {
            uint8_t byte = x & 0x7F;
            x >>= 7;
            if (x) byte |= 0x80;
            out.push_back(char(byte));
        } while (x);
    }

    static uint64_t saved_mantissa_bits(const gapped_section& sect, const gapped_recent_dists& incoming)
    {
        uint64_t saved = 0;

        for (uint8_t k = 0; k < sect.num_new; k++) {
            const uint64_t d = sect.new_dists[k];
            const auto known_end = sect.new_dists.begin() + k;
            uint64_t idx = k;

            for (uint64_t e : incoming.dist) {
                if (e == 0) break;
                if (std::find(sect.new_dists.begin(), known_end, e) != known_end) continue;
                if (idx == 3) break;

                if (e == d) {
                    saved += std::bit_width(d) - 1;
                    break;
                }

                idx++;
            }
        }

        return saved;
    }

    template <typename fnc_t>
    void for_each_char(uint64_t beg, uint64_t len, fnc_t& fnc) const
    {
        if constexpr (lce::text::is_direct_text_v<text_t>) {
            const uint8_t* data = reinterpret_cast<const uint8_t*>(T.data()) + beg;
            for (uint64_t k = 0; k < len; k++) fnc(data[k]);
        } else {
            auto cursor = T.cursor_at(beg);
            for (uint64_t k = 0; k < len; k++) fnc(T.to_char(cursor.next()));
        }
    }

    template <bool write>
    void encode_section(uint64_t s, uint16_t t)
    {
        gapped_section& sect = sections[s];
        const bool pack = write ? two_bit : scan;
        gapped_byte_writer lits(writer, stream_beg[gapped_lits] + sect.lit_beg, t, lit_write_buffer);
        gapped_bit_writer packed(writer, stream_beg[gapped_lits], 2 * sect.lit_beg, t, lit_write_buffer / 4, sect.patches);
        gapped_byte_writer excs(writer, stream_beg[gapped_excs] + sect.exc_off, t, token_write_buffer);
        gapped_byte_writer lit_lens(writer, stream_beg[gapped_lit_lens] + sect.lit_len_off, t, token_write_buffer);
        gapped_byte_writer copy_lens(writer, stream_beg[gapped_copy_lens] + sect.copy_len_off, t, token_write_buffer);
        gapped_byte_writer classes(writer, stream_beg[gapped_classes] + sect.class_off, t, token_write_buffer);
        gapped_bit_writer mantissas(writer, stream_beg[gapped_mantissas], sect.mantissa_off, t, token_write_buffer, sect.patches);
        gapped_recent_dists local;
        gapped_recent_dists global = write ? sect.incoming : gapped_recent_dists();
        uint64_t i = gapped->section_begin(s);
        uint64_t pending = 0;
        uint64_t copies = 0;
        uint64_t lit_idx = write ? sect.lit_beg : 0;
        uint64_t excs_end = write ? sect.prev_excs_end : 0;
        uint64_t run_beg = 0;
        uint64_t run_len = 0;
        uint8_t run_char = 0;

        auto end_run = [&]() {
            if (run_len == 0) return;

            if constexpr (write) {
                excs.put_vbyte(run_beg - excs_end);
                excs.put_vbyte(run_len);
                excs.put(run_char);
            } else {
                if (sect.has_excs) {
                    sect.exc_bytes += vbyte_size(run_beg - excs_end);
                } else {
                    sect.has_excs = true;
                    sect.first_exc = run_beg;
                }

                sect.exc_bytes += vbyte_size(run_len) + 1;
                sect.excs_end = run_beg + run_len;
            }

            excs_end = run_beg + run_len;
            run_len = 0;
        };

        auto add_char = [&](uint8_t c) {
            if (regular[c]) {
                end_run();
                if constexpr (!write) sect.regular++;
            } else if (run_len > 0 && run_char == c) {
                run_len++;
            } else {
                end_run();
                run_beg = lit_idx;
                run_len = 1;
                run_char = c;
            }

            if constexpr (write) packed.put(code[c], 2);
            lit_idx++;
        };

        auto add_literals = [&](uint64_t beg, uint64_t len) {
            pending += len;
            if constexpr (!write) sect.literals += len;

            if (pack) {
                for_each_char(beg, len, add_char);
            } else if constexpr (write) {
                if constexpr (lce::text::is_direct_text_v<text_t>) {
                    lits.put_bytes(T.data() + beg, len);
                } else {
                    auto cursor = T.cursor_at(beg);
                    lits.put_chars(len, [&]() { return T.to_char(cursor.next()); });
                }
            }
        };

        auto add_copy = [&](uint64_t src, uint64_t len) {
            uint64_t d = i - src;

            if (local.find(d) != 0) {
                for (uint64_t c : local.dist) {
                    if (c == 0) break;

                    if (c == d || (c <= i && T.equal(i - c, i, len))) {
                        d = c;
                        break;
                    }
                }
            }

            const int k = global.find(d);
            const uint8_t width = uint8_t(std::bit_width(d));

            if constexpr (write) {
                lit_lens.put_vbyte(copies == 0 ? sect.first_lit_len : pending);
                copy_lens.put_vbyte(len);
                classes.put(k < 0 ? uint8_t(3 + width) : uint8_t(k));
                if (k < 0) mantissas.put(d, width - 1);
            } else {
                if (copies == 0) {
                    sect.lead = pending;
                } else {
                    sect.lit_len_bytes += vbyte_size(pending);
                }

                sect.copy_len_bytes += vbyte_size(len);
                if (k < 0) sect.mantissa_bits += width - 1;
                if (local.find(d) < 0 && sect.num_new < 3) sect.new_dists[sect.num_new++] = d;
            }

            local.use(d);
            global.use(d);
            pending = 0;
            copies++;
        };

        auto add = [&](lz77_sss::factor f) {
            const uint64_t src = f.src;
            const uint64_t len = f.len;

            if (f.is_gap()) {
                add_literals(i, src);
                i += src;
            } else if (len < min_lpf_len) {
                add_literals(i, len);
                i += len;
            } else {
                add_copy(src, len);
                i += len;
            }
        };

        lz77_sss::factor_sink sink(add);
        gapped->emit_section(s, sink);
        sink.flush();
        end_run();

        if constexpr (write) {
            packed.finish();
            lits.flush();
            excs.flush();
            lit_lens.flush();
            copy_lens.flush();
            classes.flush();
            mantissas.finish();
        } else {
            sect.copies = copies;
            sect.trail = pending;
            sect.local = local;
        }
    }

    void arrange()
    {
        uint64_t lits = 0;
        uint64_t regular_lits = 0;
        uint64_t carry = 0;
        uint64_t lit_len_bytes = 0;
        uint64_t copy_len_bytes = 0;
        uint64_t copies = 0;
        uint64_t mantissa_bits = 0;
        gapped_recent_dists incoming;

        for (gapped_section& sect : sections) {
            sect.lit_beg = lits;
            sect.incoming = incoming;
            sect.copy_len_off = copy_len_bytes;
            sect.class_off = copies;
            sect.mantissa_off = mantissa_bits;
            lits += sect.literals;
            regular_lits += sect.regular;
            copy_len_bytes += sect.copy_len_bytes;
            copies += sect.copies;
            mantissa_bits += sect.mantissa_bits - saved_mantissa_bits(sect, incoming);

            if (sect.copies > 0) {
                sect.first_lit_len = carry + sect.lead;
                sect.lit_len_off = lit_len_bytes;
                lit_len_bytes += vbyte_size(sect.first_lit_len) + sect.lit_len_bytes;
                carry = sect.trail;
            } else {
                carry += sect.literals;
            }

            for (int k = 2; k >= 0; k--) {
                if (sect.local.dist[k] != 0) incoming.use(sect.local.dist[k]);
            }
        }

        two_bit = scan && lits > 0 && regular_lits >= min_two_bit_share * lits;
        uint64_t exc_bytes = 0;

        if (two_bit) {
            uint64_t excs_end = 0;

            for (gapped_section& sect : sections) {
                sect.exc_off = exc_bytes;
                sect.prev_excs_end = excs_end;
                if (!sect.has_excs) continue;
                exc_bytes += vbyte_size(sect.lit_beg + sect.first_exc - excs_end) + sect.exc_bytes;
                excs_end = sect.lit_beg + sect.excs_end;
            }
        }

        total_literals = lits;
        final_lit_len = carry;
        std::array<uint64_t, num_gapped_streams> sizes;
        sizes[gapped_lits] = two_bit ? div_ceil<uint64_t>(2 * lits, 8) : lits;
        sizes[gapped_excs] = exc_bytes;
        sizes[gapped_lit_lens] = lit_len_bytes + vbyte_size(carry);
        sizes[gapped_copy_lens] = copy_len_bytes;
        sizes[gapped_classes] = copies;
        sizes[gapped_mantissas] = div_ceil<uint64_t>(mantissa_bits, 8);
        header.clear();

        if (two_bit) {
            for (uint8_t c : chars) header.push_back(char(c));
        }

        for (uint64_t size : sizes) append_vbyte(header, size);
        stream_beg[gapped_lits] = base + header.size();
        for (uint64_t k = 1; k < num_gapped_streams; k++) stream_beg[k] = stream_beg[k - 1] + sizes[k - 1];
        final_lit_len_off = stream_beg[gapped_lit_lens] + lit_len_bytes;
    }

    const text_t& T;
    positional_writer& writer;
    uint64_t base;
    bool fasta_active;
    bool scan = false;
    bool two_bit = false;
    std::array<uint8_t, 4> chars { };
    std::array<uint8_t, 256> code { };
    std::array<bool, 256> regular { };
    const lz77_sss::gapped_factorization* gapped = nullptr;
    std::vector<gapped_section> sections;
    std::string header;
    std::array<uint64_t, num_gapped_streams> stream_beg { };
    uint64_t total_literals = 0;
    uint64_t final_lit_len = 0;
    uint64_t final_lit_len_off = 0;
};

template <typename text_t>
void encode_gapped(const text_t& T, const fasta_headers& headers, const char_histogram& histogram)
{
    uint64_t base = 0;

    {
        std::ofstream tmp_ofile(tmp_file_path, std::ios::binary);
        const uint8_t flags = 0;
        tmp_ofile.write((char*) &flags, 1);
        tmp_ofile.write((char*) &bytes_input, 8);

        if (headers.active) {
            encode_vbyte<uint64_t>(tmp_ofile, headers.size());

            for (uint64_t k = 0; k < headers.size(); k++) {
                const std::string_view line = headers.line(k);
                encode_vbyte<uint64_t>(tmp_ofile, headers.seq_pos[k]);
                encode_vbyte<uint64_t>(tmp_ofile, line.size());
                tmp_ofile.write(line.data(), line.size());
            }

            encode_vbyte<uint64_t>(tmp_ofile, T.size());
        }

        base = uint64_t(tmp_ofile.tellp());
    }

    positional_writer writer(tmp_file_path, num_threads);
    gapped_encoder<text_t> encoder(T, histogram, writer, base, headers.active);
    bool encoded = false;

    lz77_sss::factorize_gapped(T, [&](const lz77_sss::gapped_factorization& gapped) {
        encoder.encode(&gapped);
        encoded = true;
    }, { .num_threads = num_threads, .log = !quiet });

    if (!encoded) encoder.encode(nullptr);
    gaps_length = encoder.num_literals();
    two_bit_literals = encoder.uses_two_bit();
}

std::string cpu_list()
{
    std::string list;

    for (uint32_t i = 0; i < num_threads; i++) {
        if (!list.empty()) list += ",";
        list += std::to_string(num_threads == omp_get_max_threads() ? i : (2 * i));
    }

    return list;
}

postcompressor_command compress_command(uint64_t bsc_block)
{
    const std::string p = std::to_string(num_threads);
    const std::string q = std::to_string(post_compression_quality);
    const std::string in = quote(tmp_file_path);
    const std::string out = quote(output_file_path);
    const std::string ultra = post_compression_quality > 19 ? " --ultra" : "";
    using enum postcompressor_kind;

    switch (spec->kind) {
        case zstd: return { "zstd -c -q -" + q + ultra + " -T" + p + " > " + out, progress_source::pipe_input };
        case xz: return { "xz -c -q -" + q + " -T" + p + " > " + out, progress_source::pipe_input };
        case lzma: return { "xz --format=lzma -c -q -" + q + " > " + out, progress_source::pipe_input };
        case gzip: return { "gzip -c -q -" + q + " > " + out, progress_source::pipe_input };
        case pigz: return { "pigz -c -q -" + q + " -p " + p + " > " + out, progress_source::pipe_input };
        case bzip2: return { "bzip2 -c -q -" + q + " > " + out, progress_source::pipe_input };
        case pbzip2: return { "pbzip2 -c -q -" + q + " -p" + p + " > " + out, progress_source::pipe_input };
        case lbzip2: return { "lbzip2 -c -q -" + q + " -n " + p + " > " + out, progress_source::pipe_input };
        case brotli: return { "brotli -c -q " + q + " > " + out, progress_source::pipe_input };
        case lz4: return { "lz4 -c -q -" + q + " > " + out, progress_source::pipe_input };
        case lzop: return { "lzop -c -q -" + q + " > " + out, progress_source::pipe_input };
        case lzip: return { "lzip -c -q -" + q + " > " + out, progress_source::pipe_input };
        case plzip: return { "plzip -c -q -" + q + " -n " + p + " > " + out, progress_source::pipe_input };
        case bzip3: return { "bzip3 -e -c -j " + p + " -b " + q + " > " + out, progress_source::pipe_input };
        case sevenzip: return { "7z a -t7z -mx=" + q + " -mmt=" + p + " -bso0 -bsp2 " + out + " " + in + " 2>&1", progress_source::tool_stderr };
        case bsc: return { "bsc e " + in + " " + out + " -b" + std::to_string(bsc_block) + " -e2t", progress_source::tool_stdout };
    }

    return { };
}

postcompressor_command decompress_command()
{
    const std::string p = std::to_string(num_threads);
    const std::string in = quote(input_file_path);
    const std::string out = quote(tmp_file_path);
    using enum postcompressor_kind;

    switch (spec->kind) {
        case zstd: return { "zstd -d -c -q > " + out, progress_source::pipe_input };
        case xz: return { "xz -d -c -q -T" + p + " > " + out, progress_source::pipe_input };
        case lzma: return { "xz --format=lzma -d -c -q > " + out, progress_source::pipe_input };
        case gzip: return { "gzip -d -c -q > " + out, progress_source::pipe_input };
        case pigz: return { "pigz -d -c -q -p " + p + " > " + out, progress_source::pipe_input };
        case bzip2: return { "bzip2 -d -c -q > " + out, progress_source::pipe_input };
        case pbzip2: return { "pbzip2 -d -c -q -p" + p + " > " + out, progress_source::pipe_input };
        case lbzip2: return { "lbzip2 -d -c -q -n " + p + " > " + out, progress_source::pipe_input };
        case brotli: return { "brotli -d -c > " + out, progress_source::pipe_input };
        case lz4: return { "lz4 -d -c -q > " + out, progress_source::pipe_input };
        case lzop: return { "lzop -d -c -q > " + out, progress_source::pipe_input };
        case lzip: return { "lzip -d -c -q > " + out, progress_source::pipe_input };
        case plzip: return { "plzip -d -c -q -n " + p + " > " + out, progress_source::pipe_input };
        case bzip3: return { "bzip3 -d -c -j " + p + " > " + out, progress_source::pipe_input };
        case sevenzip: return { "7z x -so -bsp2 -bso0 " + in + " 2>&1 > " + out, progress_source::tool_stderr };
        case bsc: return { "bsc d " + in + " " + out, progress_source::tool_stdout };
    }

    return { };
}

void run_postcompressor(const postcompressor_command& pc, const std::string& feed_path, bool pin_cpus, const std::string& phase)
{
    std::string command;
    if (std::filesystem::exists("/usr/bin/time")) command += "/usr/bin/time -v -o " + quote(log_file_path) + " ";
    if (pin_cpus && binary_installed("taskset")) command += "taskset -c " + cpu_list() + " ";
    command += pc.command;
    if (pc.progress != progress_source::tool_stderr) command += " 2> " + quote(err_file_path);
    const bool feed = pc.progress == progress_source::pipe_input;
    log_phase_begin(!quiet, phase);
    const auto t = now();
    int status = -1;

    {
        FILE* pipe = open_pipe(command, feed);

        if (pipe != nullptr && feed) {
            #ifndef _WIN32
            auto sigpipe = std::signal(SIGPIPE, SIG_IGN);
            #endif
            const uint64_t total = std::filesystem::file_size(feed_path);
            phase_progress progress(!quiet, total);
            std::ifstream in(feed_path, std::ios::binary);
            std::string buffer;
            no_init_resize(buffer, feed_block);
            uint64_t done = 0;

            while (done < total && in.good()) {
                const uint64_t take = std::min<uint64_t>(feed_block, total - done);
                in.read(buffer.data(), take);
                if (std::fwrite(buffer.data(), 1, take, pipe) != take) break;
                done += take;
                progress.reached(done);
            }

            status = close_pipe(pipe);
            #ifndef _WIN32
            std::signal(SIGPIPE, sigpipe);
            #endif
        } else if (pipe != nullptr) {
            phase_progress progress(!quiet, 100);
            uint64_t number = 0;
            bool digits = false;

            for (int c = std::fgetc(pipe); c != EOF; c = std::fgetc(pipe)) {
                if (c >= '0' && c <= '9') {
                    number = number * 10 + uint64_t(c - '0');
                    digits = true;
                } else {
                    if (c == '%' && digits && number <= 100) progress.reached(number);
                    number = 0;
                    digits = false;
                }
            }

            status = close_pipe(pipe);
        }
    }

    if (status != 0) {
        if (!quiet) {
            std::cout << std::endl << "error: " << spec->name << " failed" << std::endl;
            std::ifstream err(err_file_path);
            std::cout << err.rdbuf() << std::flush;
        }

        std::filesystem::remove(tmp_file_path);
        std::filesystem::remove(log_file_path);
        std::filesystem::remove(err_file_path);
        exit(-1);
    }

    if (!quiet) log_runtime(t);
}

uint64_t child_peak_rss_kib()
{
    if (!std::filesystem::exists(log_file_path)) return 0;
    std::ifstream log_file(log_file_path);
    std::string log_file_str;
    uint64_t log_file_length = std::filesystem::file_size(log_file_path);
    no_init_resize(log_file_str, log_file_length);
    log_file.read(log_file_str.data(), log_file_length);
    log_file.close();
    std::string str_to_find = "Maximum resident set size (kbytes): ";
    if (log_file_str.find(str_to_find) == std::string::npos) return 0;
    uint64_t beg = log_file_str.find(str_to_find) + str_to_find.length();
    uint64_t len = log_file_str.find("\n", beg) - beg;
    return atol(log_file_str.substr(beg, len).c_str());
}

void compress()
{
    time_start = now();

    input_file.open(input_file_path);
    if (!input_file.good()) help("error: could not read <input_file>");

    if (output_file_path == "") output_file_path = input_file_path;
    output_file_path += ".ssszip." + postcompressor;
    if (std::filesystem::weakly_canonical(output_file_path) == std::filesystem::weakly_canonical(input_file_path))
        help("error: the output file must differ from <input_file>");

    bytes_input = std::filesystem::file_size(input_file_path);
    input_file.close();
    fasta_headers headers;
    char_histogram histogram { };
    uint64_t time_read = 0;

    with_text_from_file(input_file_path, bytes_input, encoding, aprx_factorization, fasta, headers,
        4 * lz77_sss::default_tau, num_threads, !quiet, [&](auto T) {
        time_read = time_diff_ns(time_start, now());

        if (result_log::path != "") {
            result_log::out.open(result_log::path, std::ios_base::app);
            result_log::text_name = input_file_path.substr(input_file_path.find_last_of("/\\") + 1);

            result_log::out << "RESULT"
                << " text_name=" << result_log::text_name
                << " type=compress"
                << " num_threads=" << num_threads
                << " n=" << bytes_input
                << " compressor=ssszip_" << postcompressor
                << " post_compression_quality=" << post_compression_quality
                << " time_read=" << time_read;
        }

        encode_gapped(T, headers, histogram);
    }, &histogram);

    double rel_len_gaps = gaps_length / (double) bytes_input;

    if (!quiet) {
        uint64_t bytes_gapped = std::filesystem::file_size(tmp_file_path);
        std::cout << "gapped factorization size = " << format_size(bytes_gapped);
        std::cout << ", relative length of the gaps = " << 100.0 * rel_len_gaps << " %";
        if (two_bit_literals) std::cout << ", 2-bit literals";
        std::cout << std::endl;
    }

    const time_point_t time_factorized = now();
    uint64_t bsc_block = uint64_t { post_compression_quality } << 20;

    if (spec->kind == postcompressor_kind::bsc && !post_compression_quality_given) {
        bsc_block = std::clamp<uint64_t>(uint64_t(malloc_count_peak() / bsc_bytes_per_block_byte),
            bsc_min_block, bsc_max_block);
    }

    const std::string setting = spec->kind == postcompressor_kind::bsc ? "block size " + format_size(bsc_block)
        : spec->kind == postcompressor_kind::bzip3 ? "block size " + std::to_string(post_compression_quality) + " MB"
        : "level " + std::to_string(post_compression_quality);
    std::filesystem::remove(output_file_path);
    run_postcompressor(compress_command(bsc_block), tmp_file_path, true,
        "compressing gapped factorization (" + std::string(spec->name) + ", " + setting + ")");
    const time_point_t time_compressed = now();

    uint64_t time_compressor = time_diff_ns(time_factorized, time_compressed);
    uint64_t gapped_peak = malloc_count_peak();
    uint64_t postcompressor_peak = child_peak_rss_kib() * 1024;
    uint64_t memory_peak = std::max(gapped_peak, postcompressor_peak);
    uint64_t time_total = time_diff_ns(time_start, time_compressed);
    bytes_compressed = std::filesystem::file_size(output_file_path);
    double compression_ratio = bytes_input / (double) bytes_compressed;
    std::filesystem::remove(log_file_path);
    std::filesystem::remove(err_file_path);
    std::filesystem::remove(tmp_file_path);

    if (!quiet) {
        std::cout << "postcompressor peak memory = " << format_size(postcompressor_peak) << std::endl;
        std::cout << "total time = " << format_time(time_total) << std::endl;
        std::cout << "total throughput = " << format_throughput(bytes_input, time_total) << std::endl;
        std::cout << "total peak memory consumption = " << format_size(memory_peak)
                  << " (" << (100.0 * memory_peak) / std::max<uint64_t>(1, bytes_input) << " % of input)" << std::endl;
        std::cout << "output file size = " << format_size(bytes_compressed) << std::endl;
        std::cout << "compression ratio = " << compression_ratio << std::endl;
    }

    if (result_log::path != "") {
        result_log::out
            << " time_compressor=" << time_compressor
            << " mem_peak_gapped=" << gapped_peak
            << " mem_peak_compressor=" << postcompressor_peak
            << " rel_len_gaps=" << rel_len_gaps
            << " two_bit_literals=" << two_bit_literals
            << " time=" << time_total
            << " throughput=" << throughput_mb_per_s(bytes_input, time_total)
            << " mem_peak=" << memory_peak
            << " bytes_comp=" << bytes_compressed
            << " comp_ratio=" << compression_ratio
            << std::endl;
    }
}

[[noreturn]] void abort_decoding(const std::string& message)
{
    if (!quiet) std::cout << std::endl << "error: " << message << std::endl;
    std::error_code ec;
    std::filesystem::remove(tmp_file_path, ec);
    std::filesystem::remove(tmp_file_path + "_seq", ec);
    std::filesystem::remove(log_file_path, ec);
    std::filesystem::remove(err_file_path, ec);
    exit(-1);
}

class gapped_byte_reader {
public:
    gapped_byte_reader(const std::string& path, uint64_t beg, uint64_t len, uint64_t capacity)
        : file(path, std::ios::binary)
        , off(beg)
        , left(len)
    {
        file.seekg(std::streamoff(beg), std::ios::beg);
        no_init_resize(buffer, std::max<uint64_t>(1, std::min(capacity, len)));
    }

    uint64_t position() const { return off - (fill - at); }

    bool done() const { return at == fill && left == 0; }

    uint8_t get()
    {
        if (at == fill) refill();
        return uint8_t(buffer[at++]);
    }

    uint64_t get_vbyte()
    {
        uint64_t x = 0;

        for (uint64_t shift = 0;; shift += 7) {
            const uint8_t byte = get();
            x |= uint64_t(byte & 0x7F) << shift;
            if ((byte & 0x80) == 0) return x;
        }
    }

    const char* next(uint64_t& len)
    {
        if (at == fill) refill();
        len = std::min(len, fill - at);
        const char* data = buffer.data() + at;
        at += len;
        return data;
    }

    void read(char* out, uint64_t len)
    {
        while (len > 0) {
            uint64_t take = len;
            const char* data = next(take);
            std::memcpy(out, data, take);
            out += take;
            len -= take;
        }
    }

private:
    void refill()
    {
        if (left == 0) abort_decoding("the gapped factorization is truncated");
        fill = std::min<uint64_t>(left, buffer.size());
        file.read(buffer.data(), std::streamsize(fill));
        if (!file) abort_decoding("could not read the gapped factorization");
        off += fill;
        left -= fill;
        at = 0;
    }

    std::ifstream file;
    std::string buffer;
    uint64_t off;
    uint64_t left;
    uint64_t fill = 0;
    uint64_t at = 0;
};

class gapped_bit_reader {
public:
    explicit gapped_bit_reader(gapped_byte_reader& in)
        : in(in)
    { }

    uint64_t get(uint8_t width)
    {
        uint64_t x = 0;

        while (width > 0) {
            if (bits == 0) {
                byte = in.get();
                bits = 8;
            }

            const uint8_t take = std::min(width, bits);
            x = (x << take) | ((byte >> (bits - take)) & ((1u << take) - 1));
            bits -= take;
            width -= take;
        }

        return x;
    }

private:
    gapped_byte_reader& in;
    uint8_t byte = 0;
    uint8_t bits = 0;
};

class gapped_two_bit_reader {
public:
    gapped_two_bit_reader(gapped_byte_reader& packed, gapped_byte_reader& excs, const std::array<char, 4>& chars)
        : packed(packed)
        , excs(excs)
        , chars(chars)
    {
        for (uint16_t b = 0; b < 256; b++) {
            for (uint8_t k = 0; k < 4; k++) quads[b][k] = chars[(b >> (6 - 2 * k)) & 3];
        }

        next_run();
    }

    void read(char* out, uint64_t len)
    {
        uint64_t k = 0;

        for (; k < len && codes > 0; k++, codes--) {
            out[k] = chars[byte >> 6];
            byte = uint8_t(byte << 2);
        }

        for (; k + 4 <= len; k += 4) std::memcpy(out + k, quads[packed.get()].data(), 4);

        if (k < len) {
            byte = packed.get();
            codes = 4;

            for (; k < len; k++, codes--) {
                out[k] = chars[byte >> 6];
                byte = uint8_t(byte << 2);
            }
        }

        const uint64_t end = index + len;

        while (run_beg < end) {
            const uint64_t from = std::max(run_beg, index);
            const uint64_t to = std::min(run_end, end);
            std::memset(out + (from - index), run_char, to - from);
            if (run_end > end) break;
            next_run();
        }

        index = end;
    }

private:
    void next_run()
    {
        if (excs.done()) {
            run_beg = std::numeric_limits<uint64_t>::max();
            return;
        }

        run_beg = run_end + excs.get_vbyte();
        run_end = run_beg + excs.get_vbyte();
        run_char = char(excs.get());
    }

    gapped_byte_reader& packed;
    gapped_byte_reader& excs;
    std::array<char, 4> chars;
    std::array<std::array<char, 4>, 256> quads;
    uint64_t index = 0;
    uint64_t run_beg = 0;
    uint64_t run_end = 0;
    uint8_t byte = 0;
    uint8_t codes = 0;
    char run_char = 0;
};

void decode_gapped()
{
    const time_point_t time_decode = now();
    const uint64_t bytes_gapped = std::filesystem::file_size(tmp_file_path);
    gapped_byte_reader header(tmp_file_path, 0, bytes_gapped, stream_read_buffer);
    const uint8_t flags = header.get();
    if ((flags & gapped_streams_flag) == 0) abort_decoding("the file has been compressed by an older version of ssszip");
    header.read((char*) &bytes_input, 8);
    const bool has_fasta = (flags & gapped_fasta_flag) != 0;
    const bool two_bit = (flags & gapped_two_bit_flag) != 0;
    fasta_headers headers;
    uint64_t bytes_decoded = bytes_input;

    if (has_fasta) {
        const uint64_t num_headers = header.get_vbyte();
        headers.seq_pos.resize(num_headers);
        headers.text_end.resize(num_headers);

        for (uint64_t k = 0; k < num_headers; k++) {
            headers.seq_pos[k] = header.get_vbyte();
            const uint64_t length = header.get_vbyte();
            const uint64_t start = headers.text.size();
            no_init_resize(headers.text, start + length);
            header.read(headers.text.data() + start, length);
            headers.text_end[k] = headers.text.size();
        }

        bytes_decoded = header.get_vbyte();
    }

    std::array<char, 4> chars { };

    if (two_bit) {
        for (char& c : chars) c = char(header.get());
    }

    std::array<uint64_t, num_gapped_streams> sizes;
    std::array<uint64_t, num_gapped_streams> begs;
    for (uint64_t& size : sizes) size = header.get_vbyte();
    begs[gapped_lits] = header.position();
    for (uint64_t k = 1; k < num_gapped_streams; k++) begs[k] = begs[k - 1] + sizes[k - 1];
    gapped_byte_reader lits(tmp_file_path, begs[gapped_lits], sizes[gapped_lits], bulk_io_size);
    gapped_byte_reader excs(tmp_file_path, begs[gapped_excs], sizes[gapped_excs], stream_read_buffer);
    gapped_byte_reader lit_lens(tmp_file_path, begs[gapped_lit_lens], sizes[gapped_lit_lens], stream_read_buffer);
    gapped_byte_reader copy_lens(tmp_file_path, begs[gapped_copy_lens], sizes[gapped_copy_lens], stream_read_buffer);
    gapped_byte_reader classes(tmp_file_path, begs[gapped_classes], sizes[gapped_classes], stream_read_buffer);
    gapped_byte_reader mantissa_bytes(tmp_file_path, begs[gapped_mantissas], sizes[gapped_mantissas], stream_read_buffer);
    gapped_bit_reader mantissas(mantissa_bytes);
    gapped_two_bit_reader unpacker(lits, excs, chars);

    const std::string sequence_file_path = has_fasta
        ? (tmp_file_path + "_seq") : output_file_path;
    file_decoder decoder(sequence_file_path, bytes_decoded, decode_in_ram, num_threads);
    if (!quiet) std::cout << "reverting gapped factorization ("
        << format_size(bytes_gapped) << ")" << std::flush;
    std::string buffer;
    if (two_bit) no_init_resize(buffer, stream_read_buffer);
    gapped_recent_dists recent;

    while (true) {
        uint64_t left = lit_lens.get_vbyte();
        if (left > bytes_decoded - decoder.position()) abort_decoding("the gapped factorization is corrupt");

        while (left > 0) {
            uint64_t take = left;

            if (two_bit) {
                take = std::min<uint64_t>(take, buffer.size());
                unpacker.read(buffer.data(), take);
                decoder.literals(buffer.data(), take);
            } else {
                const char* data = lits.next(take);
                decoder.literals(data, take);
            }

            left -= take;
        }

        if (decoder.position() == bytes_decoded) break;
        const uint64_t len = copy_lens.get_vbyte();
        const uint8_t cls = classes.get();
        uint64_t d = cls < 3 ? recent.dist[cls] : 0;
        if (cls > 3 && cls < 68) d = (uint64_t { 1 } << (cls - 4)) | mantissas.get(cls - 4);

        if (d == 0 || d > decoder.position() || len == 0 || len > bytes_decoded - decoder.position()) {
            abort_decoding("the gapped factorization is corrupt");
        }

        recent.use(d);
        decoder.copy(decoder.position() - d, len);
    }

    decoder.finish();

    if (has_fasta) {
        std::ifstream stripped_file(sequence_file_path, std::ios::binary);
        std::ofstream final_file(output_file_path, std::ios::binary);
        std::string block;
        no_init_resize(block, std::min<uint64_t>(
            std::max<uint64_t>(bytes_decoded, 1), 1024 * 1024));
        uint64_t at = 0;

        auto copy_sequence = [&](uint64_t until) {
            while (at < until) {
                const uint64_t take = std::min<uint64_t>(block.size(), until - at);
                stripped_file.read(block.data(), take);
                final_file.write(block.data(), take);
                at += take;
            }
        };

        for (uint64_t k = 0; k < headers.size(); k++) {
            const std::string_view line = headers.line(k);
            copy_sequence(headers.seq_pos[k]);
            final_file.write(line.data(), line.size());
        }

        copy_sequence(bytes_decoded);
        final_file.close();
        stripped_file.close();
        std::filesystem::remove(sequence_file_path);
    }

    uint64_t time_total = time_diff_ns(time_start, now());
    uint64_t postcompressor_peak = child_peak_rss_kib() * 1024;
    uint64_t memory_peak = std::max(postcompressor_peak, malloc_count_peak());
    double compression_ratio = bytes_input / (double) bytes_compressed;

    if (!quiet) {
        log_runtime(time_decode);
        std::cout << "total time = " << format_time(time_total) << std::endl;
        std::cout << "total throughput = " << format_throughput(bytes_input, time_total) << std::endl;
        std::cout << "total peak memory consumption = " << format_size(memory_peak)
                  << " (" << (100.0 * memory_peak) / std::max<uint64_t>(1, bytes_input) << " % of output)" << std::endl;
        std::cout << "output file size = " << format_size(bytes_input) << std::endl;
        std::cout << "compression ratio = " << compression_ratio << std::endl;
    }

    if (result_log::path != "") {
        result_log::out.open(result_log::path, std::ios_base::app);
        result_log::text_name = output_file_path.substr(output_file_path.find_last_of("/\\") + 1);

        result_log::out << "RESULT"
            << " text_name=" << result_log::text_name
            << " type=decompress"
            << " num_threads=" << num_threads
            << " n=" << bytes_input
            << " compressor=ssszip_" << postcompressor;

        if (post_compression_quality_given) result_log::out << " post_compression_quality=" << post_compression_quality;

        result_log::out
            << " time=" << time_total
            << " throughput=" << throughput_mb_per_s(bytes_input, time_total)
            << " mem_peak=" << memory_peak
            << " bytes_comp=" << bytes_compressed
            << " comp_ratio=" << compression_ratio
            << std::endl;
    }
}

void decompress()
{
    time_start = now();
    bytes_compressed = std::filesystem::file_size(input_file_path);
    if (output_file_path == "") output_file_path = input_file_path.substr(
        0, input_file_path.length() - postcompressor.length() - 8);
    if (std::filesystem::weakly_canonical(output_file_path) == std::filesystem::weakly_canonical(input_file_path))
        help("error: the output file must differ from <input_file>");
    run_postcompressor(decompress_command(), input_file_path, false,
        "decompressing gapped factorization (" + std::string(spec->name) + ", " + format_size(bytes_compressed) + ")");
    decode_gapped();
    std::filesystem::remove(tmp_file_path);
    std::filesystem::remove(log_file_path);
    std::filesystem::remove(err_file_path);
}

int main(int argc, char** argv)
{
    if (argc == 1) help("");
    num_threads = omp_get_max_threads();
    while (arg_idx < argc - 1) parse_arg(argc, argv);
    input_file_path = argv[arg_idx];
    if (input_file_path == "-h") help("");
    if (!std::filesystem::exists(input_file_path)) help("error: <input_file> does not exist");

    if (decompress_mode) {
        postcompressor = input_file_path.substr(input_file_path.find_last_of(".") + 1);
        if (input_file_path.length() < postcompressor.length() + 8 || input_file_path.substr(
            input_file_path.length() - postcompressor.length() - 8, 7) != ".ssszip")
            help("error: <input_file> does not have extension .ssszip.<postcompressor>");
    }

    spec = find_postcompressor(postcompressor);
    if (spec == nullptr) help("error: unsupported postcompressor '" + postcompressor + "'");

    if (!binary_installed(spec->binary)) {
        if (!quiet) std::cout << "error: postcompressor '" << postcompressor << "' is not installed ('" << spec->binary
            << "' not found in PATH), choose another one with -pc" << std::endl;
        exit(-1);
    }

    if (!decompress_mode && post_compression_quality_given &&
        (post_compression_quality < spec->min_quality || post_compression_quality > spec->max_quality)) {
        help("error: post-compression quality for " + postcompressor + " must be between "
            + std::to_string(spec->min_quality) + " and " + std::to_string(spec->max_quality));
    }

    if (!post_compression_quality_given) post_compression_quality = spec->default_quality;
    tmp_file_path = std::filesystem::temp_directory_path().string() + "/tmp_" + random_alphanumeric_string(10);
    log_file_path = std::filesystem::temp_directory_path().string() + "/log_" + random_alphanumeric_string(10);
    err_file_path = std::filesystem::temp_directory_path().string() + "/err_" + random_alphanumeric_string(10);
    if (decompress_mode) decompress(); else compress();
    return 0;
}