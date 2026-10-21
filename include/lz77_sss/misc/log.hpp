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

// starts the line of a phase, whose progress phase_progress shows; log_runtime ends it
inline void log_phase_begin(bool log, const std::string& message)
{
    if (!log) return;
    running_phase = message;
    phase_shows_progress = false;
    begin_progress(1);
    log_out() << message << std::flush;
}

inline void record_phase_time(const std::string& phase, uint64_t ns)
{
    if (result_log::out.is_open()) {
        result_log::out << " " << phase << "=" << ns;
    }
}
