/**
 * part of LukasNalbach/lz77-sss
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 * Copyright (c) Patrick Dinklage
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

#include <algorithm>
#include <cstdint>
#include <functional>
#include <ostream>

namespace lz77 {

struct factor {
    uintmax_t src;
    uintmax_t len;

    factor() : src(0), len(0) { }
    factor(factor&&) = default;
    factor& operator=(factor&&) = default;
    factor(const factor&) = default;
    factor& operator=(const factor&) = default;

    inline factor(char c) : src(c), len(0) { }

    inline factor(uintmax_t src, uintmax_t len) : src(src), len(len) { }

    bool operator==(const factor&) const = default;
    bool operator!=(const factor&) const = default;

    inline bool is_literal() const { return len == 0; }

    inline bool is_reference() const { return !is_literal(); }

    inline auto literal() const { return src; }

    inline size_t num_literals() const { return std::max(len, uintmax_t(1)); }

    friend std::ostream& operator<<(std::ostream& out, const factor& f)
    {
        out.write((char*) &f, sizeof(factor));
        return out;
    }
};

using emit_function = std::function<void(factor)>;

}
