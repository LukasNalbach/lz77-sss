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
#include <bit>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <string>

#include <lz77_sss/misc/log.hpp>
#include <lz77_sss/misc/utils.hpp>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

inline constexpr uint64_t progress_full = 1000;
inline constexpr uint64_t min_progress_ns = 1000000000;

inline std::string running_phase;
inline std::mutex progress_lock;
inline uint64_t progress_permille = 0;
inline uint64_t progress_drawn = 0;
inline uint64_t progress_drawn_width = 0;
inline bool phase_shows_progress = false;
inline std::chrono::steady_clock::time_point progress_started;

inline std::string format_coarse_time(uint64_t ns)
{
    if (ns < 1000000000) return "< 1 s";
    if (ns <= 10000000000) return std::to_string(ns / 1000000000) + " s";
    return format_time(ns);
}

inline bool stdout_is_terminal()
{
#ifdef _WIN32
    static const bool terminal = _isatty(_fileno(stdout)) != 0;
#else
    static const bool terminal = isatty(fileno(stdout)) != 0;
#endif
    return terminal;
}

inline void begin_progress()
{
    std::lock_guard<std::mutex> guard(progress_lock);
    progress_permille = 0;
    progress_drawn = 0;
    progress_drawn_width = 0;
    phase_shows_progress = false;
    progress_started = std::chrono::steady_clock::now();
}

inline void log_phase_begin(bool log, const std::string& message)
{
    if (!log) return;
    running_phase = message;
    begin_progress();
    std::cout << message << std::flush;
}

inline void log_progress(bool log, uint64_t permille)
{
    if (!log) return;
    std::lock_guard<std::mutex> guard(progress_lock);
    progress_permille = std::min(permille, progress_full);
    const uint64_t done = progress_permille;
    if (done == progress_drawn && phase_shows_progress) return;
    const uint64_t elapsed = time_diff_ns(progress_started);
    if (elapsed < min_progress_ns) return;
    std::string state = ' ' + std::to_string(done / 10) + '%';

    if (done >= 10 && done < progress_full) {
        state += ", elapsed: " + format_coarse_time(elapsed)
            + ", ETA " + format_coarse_time((elapsed * (progress_full - done)) / done);
    }

    if (stdout_is_terminal()) {
        const uint64_t before_width = progress_drawn_width;
        progress_drawn = done;
        phase_shows_progress = true;
        progress_drawn_width = running_phase.size() + state.size();
        const std::string blanks = progress_drawn_width < before_width
            ? std::string(before_width - progress_drawn_width, ' ') : std::string();
        std::cout << '\r' << running_phase << state << blanks << std::flush;
        return;
    }

    const uint64_t before = progress_drawn;
    progress_drawn = done;

    if (done / 100 > before / 100) {
        phase_shows_progress = true;
        std::cout << ' ' << (done / 100) * 10 << '%'
                  << (done >= 10 && done < progress_full
                         ? " (elapsed: " + format_coarse_time(elapsed)
                             + ", ETA " + format_coarse_time((elapsed * (progress_full - done)) / done) + ")"
                         : std::string())
                  << std::flush;
    }
}

inline void log_progress_done(bool log)
{
    std::lock_guard<std::mutex> guard(progress_lock);
    if (!log || !phase_shows_progress) return;
    progress_permille = progress_full;
    phase_shows_progress = false;
    if (!stdout_is_terminal()) return;
    const uint64_t width = progress_drawn_width > running_phase.size()
        ? progress_drawn_width - running_phase.size() : 0;
    progress_drawn_width = 0;
    std::cout << '\r' << running_phase << std::string(width, ' ')
              << '\r' << running_phase << std::flush;
}

class phase_progress {
public:
    phase_progress(bool log, uint64_t total, uint64_t reports = 4096)
        : active(log && total != 0), total(std::max<uint64_t>(total, 1))
    {
        const uint64_t grain = std::max<uint64_t>(1, this->total / std::max<uint64_t>(reports, 1));
        mask = (uint64_t(1) << (std::bit_width(grain) - 1)) - 1;
    }

    phase_progress(const phase_progress&) = delete;
    phase_progress& operator=(const phase_progress&) = delete;
    ~phase_progress() { log_progress_done(active); }

    void advance(uint64_t steps = 1)
    {
        if (active && steps != 0) tick(steps);
    }

    void reached(uint64_t value)
    {
        if (!active) return;
        uint64_t before = done.load(std::memory_order_relaxed);

        while (value > before
            && !done.compare_exchange_weak(before, value, std::memory_order_relaxed)) { }

        if (value > before) tick(0);
    }

    uint64_t grain() const { return mask + 1; }

private:
    void tick(uint64_t steps)
    {
        const uint64_t reached = done.fetch_add(steps, std::memory_order_relaxed) + steps;
        const uint64_t permille = std::min(progress_full, (reached * progress_full) / total);
        uint64_t before = drawn.load(std::memory_order_relaxed);

        while (permille > before) {
            if (drawn.compare_exchange_weak(before, permille, std::memory_order_relaxed)) {
                log_progress(true, permille);
                return;
            }
        }
    }

    bool active;
    uint64_t total;
    uint64_t mask;
    std::atomic<uint64_t> done = 0;
    std::atomic<uint64_t> drawn = 0;
};
