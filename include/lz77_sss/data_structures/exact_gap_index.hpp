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

#include <omp.h>

#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <libsais.h>
#include <libsais16.h>
#include <libsais40.h>
#include <util/threads.hpp>

#include <lz77_sss/data_structures/min_tree.hpp>
#include <lz77_sss/misc/utils.hpp>

template <typename text_t, typename sa_t>
class exact_gap_index {
public:
    struct occ {
        uint64_t src;
        uint64_t len;
    };

    struct no_step {
        void operator()(uint64_t) const { }
    };

    template <typename step_t = no_step>
    exact_gap_index(const text_t& text, const uint8_t* bytes, uint64_t size, uint16_t p, step_t step = { })
        : text(&text)
        , size(size)
        , sa(build_sa(bytes, size, p))
        , isa((step(sa.size() * sizeof(sa_t)), build_isa(sa, p)))
        , smaller(sa, size, p)
    {
        step(size_in_bytes() - sa.size() * sizeof(sa_t));
    }

    template <typename step_t = no_step>
        requires(std::is_same_v<sa_t, int32_t>)
    exact_gap_index(const text_t& text, const uint16_t* symbols, uint64_t size, uint16_t p, step_t step = { })
        : text(&text)
        , size(size)
        , sa(build_sa(symbols, size, p))
        , isa((step(sa.size() * sizeof(sa_t)), build_isa(sa, p)))
        , smaller(sa, size, p)
    {
        step(size_in_bytes() - sa.size() * sizeof(sa_t));
    }

    template <typename step_t = no_step>
    exact_gap_index(const text_t& text, std::vector<uint8_t>&& bytes, uint64_t size, uint16_t p, step_t step = { })
        : text(&text)
        , size(size)
        , sa(build_sa(std::move(bytes), size, p))
        , isa((step(sa.size() * sizeof(sa_t)), build_isa(sa, p)))
        , smaller(sa, size, p)
    {
        step(size_in_bytes() - sa.size() * sizeof(sa_t));
    }

    template <typename sym_t, typename step_t = no_step>
        requires(std::is_same_v<sym_t, int32_t> || std::is_same_v<sym_t, sa_int40_t>)
    exact_gap_index(const text_t& text, sym_t* symbols, uint64_t size, uint64_t sigma, uint16_t p, step_t step = { })
        : text(&text)
        , size(size)
        , sa(build_sa(symbols, size, sigma, p))
        , isa((step(sa.size() * sizeof(sa_t)), build_isa(sa, p)))
        , smaller(sa, size, p)
    {
        step(size_in_bytes() - sa.size() * sizeof(sa_t));
    }

    template <typename sym_t, typename step_t = no_step>
        requires(std::is_same_v<sym_t, int32_t> || std::is_same_v<sym_t, sa_int40_t>)
    exact_gap_index(const text_t& text, std::vector<sym_t>&& symbols, uint64_t size, uint64_t sigma, uint16_t p,
        step_t step = { })
        : text(&text)
        , size(size)
        , sa(build_sa(std::move(symbols), size, sigma, p))
        , isa((step(sa.size() * sizeof(sa_t)), build_isa(sa, p)))
        , smaller(sa, size, p)
    {
        step(size_in_bytes() - sa.size() * sizeof(sa_t));
    }

    exact_gap_index(const exact_gap_index&) = delete;
    exact_gap_index& operator=(const exact_gap_index&) = delete;

    occ lpf(uint64_t i, uint64_t limit) const
    {
        const uint64_t x = uint64_t(int64_t(isa[i]));
        occ best { 0, 0 };
        consider(smaller.previous_smaller(x), i, limit, best);
        consider(smaller.next_smaller(x), i, limit, best);
        return best;
    }

    uint64_t size_in_bytes() const
    {
        return 2 * size * sizeof(sa_t) + smaller.size_in_bytes();
    }

private:
    static std::vector<sa_t> build_sa(const uint8_t* bytes, uint64_t size, uint16_t p)
    {
        std::vector<sa_t> sa;
        no_init_resize(sa, size);

        const int threads = lce::util::sais_threads(p);

        if constexpr (std::is_same_v<sa_t, int32_t>) {
            if (libsais_omp(bytes, sa.data(), int32_t(size), 0, nullptr, threads) != 0) {
                throw std::runtime_error("libsais_omp failed");
            }
        } else {
            if (libsais40_impl::libsais40_omp(bytes, sa.data(), int64_t(size), 0, nullptr, threads) != 0) {
                throw std::runtime_error("libsais40_omp failed");
            }
        }

        return sa;
    }

    static std::vector<sa_t> build_sa(const uint16_t* symbols, uint64_t size, uint16_t p)
    {
        std::vector<sa_t> sa;
        no_init_resize(sa, size);

        if (libsais16_omp(symbols, sa.data(), int32_t(size), 0, nullptr, lce::util::sais_threads(p)) != 0) {
            throw std::runtime_error("libsais16_omp failed");
        }

        return sa;
    }

    static std::vector<sa_t> build_sa(std::vector<uint8_t>&& bytes, uint64_t size, uint16_t p)
    {
        std::vector<uint8_t> owned = std::move(bytes);
        return build_sa(owned.data(), size, p);
    }

    template <typename sym_t>
    static std::vector<sa_t> build_sa(sym_t* symbols, uint64_t size, uint64_t sigma, uint16_t p)
    {
        std::vector<sa_t> sa;
        no_init_resize(sa, size);

        const int threads = lce::util::sais_threads(p);

        if constexpr (std::is_same_v<sa_t, int32_t>) {
            static_assert(std::is_same_v<sym_t, int32_t>);

            if (libsais_int_omp(symbols, sa.data(), int32_t(size), int32_t(std::max<uint64_t>(sigma, 1)), 0, threads) != 0) {
                throw std::runtime_error("libsais_int_omp failed");
            }
        } else {
            static_assert(std::is_same_v<sym_t, sa_int40_t>);

            if (libsais40_impl::libsais40_long_omp(symbols, sa.data(), int64_t(size),
                    int64_t(std::max<uint64_t>(sigma, 1)), 0, threads) != 0) {
                throw std::runtime_error("libsais40_long_omp failed");
            }
        }

        return sa;
    }

    template <typename sym_t>
    static std::vector<sa_t> build_sa(std::vector<sym_t>&& symbols, uint64_t size, uint64_t sigma, uint16_t p)
    {
        std::vector<sym_t> owned = std::move(symbols);
        return build_sa(owned.data(), size, sigma, p);
    }

    static std::vector<sa_t> build_isa(const std::vector<sa_t>& sa, uint16_t p)
    {
        std::vector<sa_t> isa;
        no_init_resize(isa, sa.size());

        #pragma omp parallel for num_threads(p) schedule(dynamic, 65536)
        for (uint64_t x = 0; x < sa.size(); x++) {
            isa[uint64_t(int64_t(sa[x]))] = sa_t(int64_t(x));
        }

        return isa;
    }

    inline void consider(uint64_t y, uint64_t i, uint64_t limit, occ& best) const
    {
        if (y == smaller.none()) return;
        const uint64_t src = uint64_t(int64_t(sa[y]));
        const uint64_t len = text->lce(src, i, limit);

        if (len > best.len || (len == best.len && len > 0 && src > best.src)) {
            best = occ { src, len };
        }
    }

    const text_t* text;
    uint64_t size;
    std::vector<sa_t> sa;
    std::vector<sa_t> isa;
    min_tree<std::vector<sa_t>> smaller;
};
