/**
 * lz77/factor.hpp
 * part of pdinklag/lz77
 * 
 * MIT License
 * 
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

#ifndef _LZ77_FACTOR_HPP
#define _LZ77_FACTOR_HPP

#include <algorithm>
#include <cstdint>

namespace lz77 {

/**
 * \brief Represents a Lempel-Ziv 77 factor
 * 
 * This structure serves as a simple contract for communicating LZ77 factors, which are commonly represented as tuples of two integers.
 * A factor describes either a reference (copy \c len characters from \c src positions ago) or a literal factor (if \c len equals zero, \c src contains a character).
 * The different meanings of \ref src and \ref len depend on the context.
 * There are several queries, such as \ref is_literal , that reflect common usecases.
 * 
 * This structure is \em not meant to store LZ77 factors in a space-efficient manner.
 */
struct Factor {
    /**
     * \brief The copy source of a referencing, or the character value of a literal factor
     */
    uintmax_t src;

    /**
     * \brief The length of the referencing factor, or zero to indicate that this is a literal factor
     */
    uintmax_t len;

    Factor() : src(0), len(0) {}
    Factor(Factor&&) = default;
    Factor& operator=(Factor&&) = default;
    Factor(Factor const&) = default;
    Factor& operator=(Factor const&) = default;

    /**
     * \brief Constructs a literal factor
     * 
     * The literal value will be stored in \ref src and \ref len will be initialized as zero.
     * 
     * \param c the literal value
     */
    inline Factor(char const c) : src(c), len(0) {}

    /**
     * \brief Constructs a factor
     * 
     * \param _src the source value
     * \param _len the length value
     */
    inline Factor(uintmax_t const src, uintmax_t const len) : src(src), len(len) {}

    bool operator==(Factor const&) const = default;
    bool operator!=(Factor const&) const = default;

    /**
     * \brief Tests whether this factor is a literal factor
     * 
     * \return true if \ref len equals zero
     * \return false otherwise
     */
    inline bool is_literal() const { return len == 0; }

    /**
     * \brief Tests whether this factor is a referencing factor
     * 
     * \return true if \ref len is larger than zero
     * \return false otherwise
     */
    inline bool is_reference() const { return !is_literal(); }

    /**
     * \brief Extracts the literal value of the factor
     * 
     * This is only meaningful if \ref is_literal reports \c true .
     * 
     * \return the literal value of the factor
     */
    inline auto literal() const { return src; }

    /**
     * \brief Reports the number of literals encoded by this factor
     * 
     * For referencing factors, this equals their length.
     * For literal factors, this equals one.
     * 
     * \return the number of literals encoded by this factor
     */
    inline size_t num_literals() const { return std::max(len, uintmax_t(1)); }

    friend std::ostream& operator<<(std::ostream& out, const Factor& f) {
        out.write((char*)&f,sizeof(Factor));
        return out;
    }
};

}

#endif
