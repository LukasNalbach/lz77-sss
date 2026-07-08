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

#include <bit>
#include <fstream>
#include <vector>
#include <lz77_sss/lz77_sss.hpp>
#include <lz77_sss/misc/huffman.hpp>

std::fstream input_file;
std::fstream output_file;

int main(int argc, char** argv)
{
    if (argc != 3 && argc != 4) {
        std::cout << "usage: lz77_sss_decode <input_file> <output_file> [throughput_threshold_MB/s]" << std::endl;
        exit(-1);
    }

    input_file.open(argv[1], std::ios::in | std::ios::binary);
    uint64_t input_file_size = std::filesystem::file_size(argv[1]);
    std::cout << "input file size: " << format_size(input_file_size) << std::endl;

    if (!input_file.good()) {
        std::cout << "error: could not read <input_file>" << std::endl;
        exit(-1);
    }

    std::string output_file_name = argv[2];
    if (std::filesystem::exists(output_file_name)) std::filesystem::remove(output_file_name);
    output_file.open(output_file_name, std::ios::in | std::ios::out | std::ios::app | std::ios::binary);

    if (!output_file.good()) {
        std::cout << "error: could not write to <output_file>" << std::endl;
        exit(-1);
    }

    uint64_t n;
    input_file.read((char*) &n, 5);
    bit_reader reader(input_file);
    huffman len_huff, dist_huff;
    using factor = lz77_sss<uint64_t>::factor;
    huff_factor_iterator<factor> it(reader, len_huff, dist_huff, n);

    factor f;
    uint64_t pos_output = 0;
    std::vector<char> ring;
    uint64_t buff_size = 0;
    uint64_t buff_mask = 0;
    uint64_t buffered_start = 0;
    const uint64_t init_buff_size = uint64_t(1) << 16;
    const uint64_t max_buff_size = std::max(init_buff_size, std::bit_floor(n));
    uint64_t file_buff_size = std::max<uint64_t>(32 * 1024, n / 1000);
    double throughput_threshold = argc == 4 ? std::stod(argv[3]) : 20.0;
    const uint64_t copy_chunk = uint64_t(1) << 20;
    std::vector<char> copy_buff;
    std::string file_buff;

    auto grow_buffer = [&]() -> bool {
        uint64_t new_size = buff_size == 0 ? init_buff_size : 2 * buff_size;
        if (new_size > max_buff_size) new_size = max_buff_size;
        if (new_size == buff_size) return false;

        uint64_t new_mask = new_size - 1;
        std::vector<char> new_ring;
        no_init_resize(new_ring, new_size);
        uint64_t valid_len = buff_size == 0 ? 0
            : std::min(buff_size, pos_output - buffered_start);
        buffered_start = pos_output - valid_len;

        for (uint64_t p = buffered_start; p < pos_output; p++)
            new_ring[p & new_mask] = ring[p & buff_mask];

        ring.swap(new_ring);
        buff_size = new_size;
        buff_mask = new_mask;
        return true;
    };

    auto emit_literal = [&](char c) {
        if (buff_size) ring[pos_output & buff_mask] = c;
        output_file.write(&c, 1);
        pos_output++;
    };

    auto copy_from_ring = [&](uint64_t from, uint64_t len) {
        for (uint64_t off = 0; off < len;) {
            uint64_t chunk = std::min(len - off, copy_chunk);
            if (copy_buff.size() < chunk) copy_buff.resize(chunk);

            for (uint64_t i = 0; i < chunk; i++) {
                char c = ring[(from + i) & buff_mask];
                ring[(pos_output + i) & buff_mask] = c;
                copy_buff[i] = c;
            }

            output_file.write(copy_buff.data(), chunk);
            from += chunk;
            pos_output += chunk;
            off += chunk;
        }
    };

    auto copy_via_file = [&](uint64_t from, uint64_t len) {
        uint64_t to = pos_output;
        uint64_t chunk_cap = std::min(file_buff_size, to - from);
        if (file_buff.size() < chunk_cap) file_buff.resize(chunk_cap);

        for (uint64_t off = 0; off < len;) {
            uint64_t chunk = std::min(len - off, chunk_cap);
            output_file.seekg(from + off, std::ios::beg);
            output_file.read(file_buff.data(), chunk);
            output_file.seekp(to + off, std::ios::beg);
            output_file.write(file_buff.data(), chunk);

            if (buff_size) {
                for (uint64_t i = 0; i < chunk; i++)
                    ring[(to + off + i) & buff_mask] = file_buff[i];
            }

            off += chunk;
        }

        pos_output += len;
    };

    std::cout << "decoding (" << format_size(n) << ")" << std::flush;
    auto t1 = now();

    const uint64_t base_interval = std::clamp<uint64_t>(n / 16, uint64_t(1) << 20, uint64_t(1) << 23);
    const uint64_t min_eval_size = uint64_t(1) << 22;
    const double min_hit_gain = 0.03;
    const uint32_t max_stale = 2;

    bool measuring = true;
    bool eval_pending = false;
    bool frozen = false;
    uint32_t stale = 0;
    double hit_before_grow = 0.0;
    uint64_t seg_start_pos = 0, seg_copy = 0, seg_ring = 0;
    auto seg_start_time = t1;
    uint64_t next_check = base_interval;

    while (pos_output < n) {
        f = *it++;

        if (f.len == 0) {
            emit_literal((char) f.src);
        } else {
            seg_copy += f.len;
            if (buff_size && f.src >= buffered_start && pos_output - f.src <= buff_size) {
                copy_from_ring(f.src, f.len);
                seg_ring += f.len;
            } else {
                copy_via_file(f.src, f.len);
            }
        }

        if (pos_output >= next_check) {
            if (!measuring) {
                measuring = true;
                seg_start_pos = pos_output;
                seg_start_time = now();
                seg_copy = seg_ring = 0;
                next_check = pos_output + base_interval;
            } else {
                uint64_t elapsed = time_diff_ns(seg_start_time, now());
                double tp = elapsed == 0 ? throughput_threshold + 1.0
                    : throughput(pos_output - seg_start_pos, elapsed);
                double hit = seg_copy == 0 ? 1.0 : seg_ring / (double) seg_copy;

                if (eval_pending) {
                    if (buff_size >= min_eval_size && hit - hit_before_grow < min_hit_gain) {
                        if (++stale >= max_stale) frozen = true;
                    } else stale = 0;
                    eval_pending = false;
                }

                if (!frozen && tp < throughput_threshold) {
                    hit_before_grow = hit;
                    if (grow_buffer()) {
                        eval_pending = true;
                        measuring = false;
                        next_check = pos_output + buff_size;
                    } else {
                        frozen = true;
                    }
                }

                if (frozen) {
                    next_check = n;
                } else if (measuring) {
                    seg_start_pos = pos_output;
                    seg_start_time = now();
                    seg_copy = seg_ring = 0;
                    next_check = pos_output + base_interval;
                }
            }
        }
    }

    auto t2 = now();
    log_runtime(t1, t2);
    std::cout << "throughput: " << format_throughput(n, time_diff_ns(t1, t2)) << std::endl;
    if (buff_size > 0) std::cout << "ring buffer size: " << format_size(buff_size) << std::endl;
    std::cout << "peak memory consumption: " << format_size(malloc_count_peak()) << std::endl;
    std::cout << "compression ratio: " << n / (double) input_file_size << std::endl;
    return 0;
}
