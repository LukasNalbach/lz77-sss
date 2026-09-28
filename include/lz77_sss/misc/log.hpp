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

#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <lz77_sss/misc/utils.hpp>

namespace result_log {

inline std::string text_name;
inline std::string path;
inline std::ofstream out;
inline bool write_rows = false;
inline std::ofstream queries_out;

}

inline std::string format_time(uint64_t ns)
{
    const uint64_t s = ns / 1000000000;

    if (ns <= 10000) return std::to_string(ns) + " ns";
    if (ns <= 10000000) return std::to_string(ns / 1000) + " us";
    if (ns <= 10000000000) return std::to_string(ns / 1000000) + " ms";
    if (s < 60) return std::to_string(s) + " s";
    if (s < 3600) return std::to_string(s / 60) + " min " + std::to_string(s % 60) + " s";
    if (s < 86400) return std::to_string(s / 3600) + " h " + std::to_string(s % 3600 / 60) + " min";
    return std::to_string(s / 86400) + " d " + std::to_string(s % 86400 / 3600) + " h";
}

inline std::string format_size(uint64_t bytes)
{
    std::string size_str;

    if (bytes > 10000000000) {
        size_str = std::to_string(bytes / 1000000000) + " GB";
    } else if (bytes > 10000000) {
        size_str = std::to_string(bytes / 1000000) + " MB";
    } else if (bytes > 10000) {
        size_str = std::to_string(bytes / 1000) + " kB";
    } else {
        size_str = std::to_string(bytes) + " B";
    }

    return size_str;
}

inline double throughput_mb_per_s(uint64_t bytes, uint64_t ns) { return 1'000.0 * (bytes / (double)ns); }

inline std::string format_throughput(uint64_t bytes, uint64_t ns)
{
    return std::to_string(throughput_mb_per_s(bytes, ns)) + " MB/s";
}

inline std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2)
{
    std::cout << ", in ~ " << format_time(time_diff_ns(t1, t2)) << std::endl;
    return std::chrono::steady_clock::now();
}

inline std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t)
{
    return log_runtime(t, std::chrono::steady_clock::now());
}

inline void record_phase_time(const std::string& phase, uint64_t ns)
{
    if (result_log::out.is_open()) {
        result_log::out << " " << phase << "=" << ns;
    }
}
