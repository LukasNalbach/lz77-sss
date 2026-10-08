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

#include <cstdint>
#include <random>
#include <vector>

#include <lz77_sss/misc/utils.hpp>
#include <text/direct_text.hpp>

template <uint8_t mers_exp = 61, typename text_t = lce::text::direct_text<char>>
class rabin_karp_substring {
protected:
    static_assert(mers_exp == 31 || mers_exp == 61);
    static_assert(mers_exp == 61 || text_t::is_byte_text);

    using fp_t = std::conditional_t<mers_exp == 31, uint32_t, uint64_t>;
    using symbol_t = typename text_t::symbol_type;
    using fp_wide_t = std::conditional_t<mers_exp == 31, uint64_t, uint128_t>;

    static constexpr fp_t mersenne_prime = (fp_t(1) << mers_exp) - 1;
    static constexpr fp_wide_t mersenne_prime_sq = fp_wide_t(mersenne_prime) * mersenne_prime;

    text_t T;
    uint64_t n;
    uint64_t sqrt_n;
    fp_wide_t base;
    uint64_t sample_rate;
    uint64_t w;
    fp_wide_t pop_prec[256] = { };
    fp_wide_t pop_mul = 0;
    std::vector<fp_t> prefix_fps;
    std::vector<fp_t> base_pow_leq_sqrt_n;
    std::vector<fp_t> base_pow_step_sqrt_n;
    
    inline static fp_t mod(const fp_wide_t val)
    {
        fp_wide_t res_wide = (val >> mers_exp) + (val & mersenne_prime);
        res_wide = (res_wide >> mers_exp) + (res_wide & mersenne_prime);
        fp_t res = (fp_t)res_wide;
        res -= mersenne_prime & -(fp_t)(res >= mersenne_prime);
        return res;
    }

    inline fp_t base_pow(const uint64_t exp) const
    {
        const uint64_t exp_sqrt_n = exp / sqrt_n;
        const uint64_t offs_exp = exp - exp_sqrt_n * sqrt_n;
        return mod(fp_wide_t(base_pow_step_sqrt_n[exp_sqrt_n]) * fp_wide_t(base_pow_leq_sqrt_n[offs_exp]));
    }

public:
    rabin_karp_substring() = default;

    rabin_karp_substring(
        const char* T,
        const uint64_t n,
        const uint64_t sample_rate = 16,
        const uint64_t w = 0,
        const uint16_t p = 1)
        requires(lce::text::is_direct_text_v<text_t> && text_t::is_byte_text)
        : rabin_karp_substring(text_t(const_cast<char*>(T), n), n, sample_rate, w, p) { }

    rabin_karp_substring(
        const text_t& T,
        const uint64_t n,
        const uint64_t sample_rate = 16,
        const uint64_t w = 0,
        const uint16_t p = 1)
        : T(T)
        , n(n), sample_rate(sample_rate), w(w)
    {
        sqrt_n = std::ceil(std::sqrt(double(n)));
        std::random_device rd;
        std::mt19937_64 mt(rd());

        base = std::uniform_int_distribution<fp_t>(257, mersenne_prime - 1)(mt);
        const uint64_t num_blks = div_ceil<uint64_t>(n, sample_rate);
        no_init_resize(base_pow_leq_sqrt_n, sqrt_n + 1);
        no_init_resize(base_pow_step_sqrt_n, sqrt_n + 1);
        base_pow_leq_sqrt_n[0] = 1;
        base_pow_step_sqrt_n[0] = 1;

        for (uint64_t i = 1; i <= sqrt_n; i++) {
            base_pow_leq_sqrt_n[i] = mod(fp_wide_t(base_pow_leq_sqrt_n[i - 1]) * base);
        }

        for (uint64_t i = 1; i <= sqrt_n; i++) {
            base_pow_step_sqrt_n[i] = mod(
                fp_wide_t(base_pow_step_sqrt_n[i - 1]) *
                fp_wide_t(base_pow_leq_sqrt_n[sqrt_n]));
        }

        if (w != 0) {
            const fp_wide_t max_exp_excl = base_pow(w);
            pop_mul = max_exp_excl;

            for (fp_wide_t c = 0; c < 256; c++) {
                pop_prec[c] = mersenne_prime_sq - max_exp_excl * c;
            }
        }

        if (sample_rate >= n) {
            no_init_resize(prefix_fps, num_blks + 1);
            prefix_fps[0] = 0;
            fp_t fp = 0;

            for (uint64_t blk = 0; blk < num_blks; blk++) {
                uint64_t beg = blk * sample_rate;
                uint64_t end = std::min<uint64_t>(beg + sample_rate, n);
                fp = concat(fp, substr_fp_naive<RIGHT>(beg, end - beg), end - beg);
                prefix_fps[blk + 1] = fp;
            }

            return;
        }

        uint64_t blks_per_thr = num_blks / p;
        std::vector<fp_t> blk_fps;
        no_init_resize(blk_fps, p + 1);
        blk_fps[0] = 0;

        #pragma omp parallel num_threads(p)
        {
            uint16_t i_p = omp_get_thread_num();
            const uint64_t beg = i_p * sample_rate * blks_per_thr;
            const uint64_t end = i_p == p - 1 ? n : ((i_p + 1) * sample_rate * blks_per_thr);
            blk_fps[i_p + 1] = substr_fp_naive(beg, end - beg);
        }

        for (uint16_t i_p = 1; i_p < p; i_p++) {
            blk_fps[i_p] = concat(blk_fps[i_p - 1], blk_fps[i_p], sample_rate * blks_per_thr);
        }
        
        no_init_resize(prefix_fps, num_blks + 1);
        prefix_fps[0] = 0;

        #pragma omp parallel num_threads(p)
        {
            uint16_t i_p = omp_get_thread_num();
            const uint64_t blk_beg = i_p * blks_per_thr;
            const uint64_t blk_end = i_p == p - 1 ? num_blks : (blk_beg + blks_per_thr);
            fp_t fp = blk_fps[i_p];
            auto cursor = T.cursor_at(blk_beg * sample_rate);

            for (uint64_t blk = blk_beg; blk < blk_end;) {
                uint64_t beg = blk * sample_rate;
                uint64_t end = beg + sample_rate;
                if (blk + 1 == num_blks) [[unlikely]] end = n;
                for (uint64_t i = beg; i < end; i++) fp = push(fp, cursor.next());
                prefix_fps[++blk] = fp;
            }
        }
    }

    inline uint64_t window_size() const { return w; }

    inline uint64_t size_in_bytes() const
    {
        uint64_t size = sizeof(*this);
        size += base_pow_leq_sqrt_n.size() * sizeof(fp_t);
        size += base_pow_step_sqrt_n.size() * sizeof(fp_t);
        size += prefix_fps.size() * sizeof(fp_t);
        return size;
    }

    inline fp_t push(const fp_t fp, const symbol_t chr) const
    {
        return mod(((base * fp_wide_t(fp)) + mersenne_prime_sq) + fp_wide_t(chr));
    }

    inline fp_t roll(const fp_t fp, const symbol_t pop, const symbol_t chr) const
    {
        if constexpr (text_t::is_byte_text) {
            return mod(((base * fp_wide_t(fp)) + pop_prec[pop]) + fp_wide_t(chr));
        } else {
            const fp_wide_t pop_red = sizeof(symbol_t) < 8 ? fp_wide_t(pop) : fp_wide_t(mod(fp_wide_t(pop)));
            return mod(((base * fp_wide_t(fp)) + (mersenne_prime_sq - pop_mul * pop_red)) + fp_wide_t(chr));
        }
    }

    inline fp_t concat(const fp_t fp_left, const fp_t fp_right, const uint64_t len_right) const
    {
        return mod(fp_wide_t(base_pow(len_right)) * fp_wide_t(fp_left) + fp_wide_t(fp_right));
    }

    template <direction dir = RIGHT>
    inline fp_t substr_fp_naive(uint64_t pos, const uint64_t len) const
    {
        if constexpr (dir == LEFT) pos -= len - 1;
        fp_t fp = 0;
        auto cursor = T.cursor_at(pos);

        for (uint64_t i = 0; i < len; i++) {
            fp = push(fp, cursor.next());
        }

        return fp;
    }

    inline fp_t prefix_fp(const uint64_t pos) const
    {
        uint64_t blk = pos / sample_rate;
        uint64_t blk_beg = blk * sample_rate;
        uint64_t offs = pos - blk_beg;
        if (offs == 0) [[unlikely]] return prefix_fps[blk];
        return concat(prefix_fps[blk], substr_fp_naive<RIGHT>(blk_beg, offs), offs);
    }

    inline fp_t substr_fp_of_prefixes(const fp_t fp_beg, const fp_t fp_end, const uint64_t len) const
    {
        const fp_t fp_beg_shft = mod(fp_wide_t(fp_beg) * fp_wide_t(base_pow(len)));
        return fp_end >= fp_beg_shft ? fp_end - fp_beg_shft : mersenne_prime - (fp_beg_shft - fp_end);
    }

    template <direction dir = RIGHT>
    inline fp_t substr_fp(uint64_t pos, const uint64_t len) const
    {
        if constexpr (dir == LEFT) pos -= len - 1;
        return substr_fp_of_prefixes(prefix_fp(pos), prefix_fp(pos + len), len);
    }

    template <direction dir = RIGHT>
    inline void substr_fps(const uint64_t pos, const uint64_t* lens, const uint64_t num, uint64_t* fps) const
    {
        if (num == 0) return;
        const uint64_t first = dir == RIGHT ? pos : pos + 1 - lens[num - 1];
        uint64_t cur = first - first % sample_rate;
        fp_t fp = prefix_fps[cur / sample_rate];
        auto cursor = T.cursor_at(cur);

        auto prefix_fp_at = [&](const uint64_t target) {
            const uint64_t blk_beg = target - target % sample_rate;

            if (cur + sample_rate < blk_beg) {
                cur = blk_beg;
                fp = prefix_fps[blk_beg / sample_rate];
                cursor = T.cursor_at(std::min<uint64_t>(blk_beg, n - 1));
            }

            for (; cur < target; cur++) fp = push(fp, cursor.next());
            return fp;
        };

        if constexpr (dir == RIGHT) {
            const fp_t fp_beg = prefix_fp_at(pos);

            for (uint64_t k = 0; k < num; k++) {
                fps[k] = substr_fp_of_prefixes(fp_beg, prefix_fp_at(pos + lens[k]), lens[k]);
            }
        } else {
            for (uint64_t k = num; k-- > 0;) fps[k] = prefix_fp_at(pos + 1 - lens[k]);
            const fp_t fp_end = prefix_fp_at(pos + 1);

            for (uint64_t k = 0; k < num; k++) {
                fps[k] = substr_fp_of_prefixes(fp_t(fps[k]), fp_end, lens[k]);
            }
        }
    }
};
