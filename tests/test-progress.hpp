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

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <string>
#include <vector>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include <omp.h>

inline uint64_t fuzz_iterations(uint64_t default_iterations = 1000)
{
    if (const char* env = std::getenv("LZ77_SSS_TEST_ITERATIONS")) {
        uint64_t value = std::strtoull(env, nullptr, 10);
        if (value > 0) return value;
    }

    return default_iterations;
}

struct fuzz_functionality {
    std::string name;
    std::function<void(uint64_t)> run;
    bool parallel = false;
};

inline void run_fuzz(
    const std::string& structure, const std::vector<fuzz_functionality>& functionalities,
    uint64_t iterations = fuzz_iterations())
{
    using clock = std::chrono::steady_clock;
#ifdef _WIN32
    const bool tty = _isatty(_fileno(stderr)) != 0;
#else
    const bool tty = isatty(fileno(stderr)) != 0;
#endif
    const char* green = tty ? "\033[0;32m" : "";
    const char* reset = tty ? "\033[m" : "";

    const uint64_t total_iterations = std::max<uint64_t>(1, iterations);
    const size_t num_functionalities = std::max<size_t>(1, functionalities.size());

    auto seconds_since = [](clock::time_point t) {
        return std::chrono::duration<double>(clock::now() - t).count();
    };

    int last_pct = -1;
    auto draw = [&](const std::string& label, uint64_t done, double elapsed) {
        constexpr int width = 10;
        const double fraction = (double) done / total_iterations;
        int filled = (int) (fraction * width + 0.5);
        if (filled > width) filled = width;
        std::string bar(filled, '=');
        bar.resize(width, ' ');
        std::fprintf(stderr, "%s%s[%s]%s %-28s %7llu/%-7llu %5.1fs   ",
            tty ? "\r" : "\n", green, bar.c_str(), reset, (structure + "/" + label).c_str(),
            (unsigned long long) done, (unsigned long long) total_iterations, elapsed);
        std::fflush(stderr);
    };

    uint64_t completed = 0;
    for (size_t f = 0; f < functionalities.size(); f++) {
        const fuzz_functionality& functionality = functionalities[f];
        const uint64_t slice_iterations = total_iterations / num_functionalities
            + (f < total_iterations % num_functionalities ? 1 : 0);
        const auto slice_start = clock::now();

        if (functionality.parallel) {
            std::atomic<uint64_t> finished { 0 };

            #pragma omp parallel for schedule(dynamic)
            for (int64_t iteration = 0; iteration < (int64_t) slice_iterations; iteration++) {
                functionality.run((uint64_t) iteration);
                const uint64_t done = completed + finished.fetch_add(1) + 1;
                if (omp_get_thread_num() == 0) {
                    const int pct = (int) (100.0 * done / total_iterations + 0.5);
                    if (pct != last_pct && (tty || pct % 10 == 0)) {
                        last_pct = pct;
                        draw(functionality.name, done, seconds_since(slice_start));
                    }
                }
            }
        } else {
            for (uint64_t iteration = 0; iteration < slice_iterations; iteration++) {
                functionality.run(iteration);
                const uint64_t done = completed + iteration + 1;
                const int pct = (int) (100.0 * done / total_iterations + 0.5);
                if (pct != last_pct && (tty || pct % 10 == 0)) {
                    last_pct = pct;
                    draw(functionality.name, done, seconds_since(slice_start));
                }
            }
        }

        completed += slice_iterations;
        draw(functionality.name, completed, seconds_since(slice_start));
        std::fputc('\n', stderr);
        last_pct = -1;
    }
}
