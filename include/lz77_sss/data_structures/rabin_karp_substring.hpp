#pragma once

#include <cstdint>
#include <random>
#include <vector>

#include <lz77_sss/misc/utils.hpp>

template <uint8_t mers_exp = 61, typename pos_t = uint32_t>
class rabin_karp_substring {
protected:
    static_assert(mers_exp == 31 || mers_exp == 61);

    using fp_t = std::conditional_t<mers_exp == 31, uint32_t, uint64_t>;
    using fp_conc_t = std::conditional_t<mers_exp == 31, uint64_t, uint128_t>;

    static constexpr fp_t mers_prim = (fp_t(1) << mers_exp) - 1;
    static constexpr fp_conc_t sqr_mers_prim = fp_conc_t(mers_prim) * mers_prim;

    const uint8_t* T = nullptr;
    pos_t n;
    pos_t sqrt_n;
    fp_conc_t b;
    pos_t s;
    pos_t w;
    fp_conc_t pop_prec[256] = { };
    std::vector<fp_t> fps;
    std::vector<fp_t> b_pow_leq_sqrt_n;
    std::vector<fp_t> b_pow_step_sqrt_n;

    inline static fp_t mod(const fp_conc_t val)
    {
        const fp_conc_t v = val + 1;
        const fp_t z = ((v >> mers_exp) + v) >> mers_exp;
        return (val + z) & mers_prim;
    }

    inline fp_t b_pow(const pos_t exp) const
    {
        const pos_t exp_sqrt_n = exp / sqrt_n;
        const pos_t offs_exp = exp - exp_sqrt_n * sqrt_n;
        return mod(fp_conc_t(b_pow_step_sqrt_n[exp_sqrt_n]) * fp_conc_t(b_pow_leq_sqrt_n[offs_exp]));
    }

public:
    rabin_karp_substring() = default;

    template <typename char_t>
    rabin_karp_substring(
        const char_t* T,
        const pos_t n,
        const pos_t s = 16,
        const pos_t w = 0,
        const uint16_t p = 1)
        : T(reinterpret_cast<const uint8_t*>(T))
        , n(n), s(s), w(w)
    {
        static_assert(sizeof(char_t) == 1);

        sqrt_n = std::ceil(std::sqrt(double(n)));
        std::random_device rd;
        std::mt19937_64 mt(rd());
        std::uniform_int_distribution<fp_t> distrib(257, mers_prim);
        b = distrib(mt);
        const pos_t num_blks = div_ceil<pos_t>(n, s);

        no_init_resize(b_pow_leq_sqrt_n, sqrt_n + 1);
        no_init_resize(b_pow_step_sqrt_n, sqrt_n + 1);
        b_pow_leq_sqrt_n[0] = 1;
        b_pow_step_sqrt_n[0] = 1;

        for (pos_t i = 1; i <= sqrt_n; i++) {
            b_pow_leq_sqrt_n[i] = mod(fp_conc_t(b_pow_leq_sqrt_n[i - 1]) * b);
        }

        for (pos_t i = 1; i <= sqrt_n; i++) {
            b_pow_step_sqrt_n[i] = mod(
                fp_conc_t(b_pow_step_sqrt_n[i - 1]) *
                fp_conc_t(b_pow_leq_sqrt_n[sqrt_n]));
        }

        if (w != 0) {
            const fp_conc_t max_exp_excl = b_pow(w);

            for (fp_conc_t c = 0; c < 256; c++) {
                pop_prec[c] = sqr_mers_prim - max_exp_excl * c;
            }
        }

        if (s >= n) return;
        pos_t blks_per_thr = num_blks / p;
        std::vector<fp_t> blk_fps;
        no_init_resize(blk_fps, p + 1);
        blk_fps[0] = 0;

        #pragma omp parallel num_threads(p)
        {
            uint16_t i_p = omp_get_thread_num();
            const pos_t beg = i_p * s * blks_per_thr;
            const pos_t end = i_p == p - 1 ? n : ((i_p + 1) * s * blks_per_thr);
            blk_fps[i_p + 1] = substring_naive(beg, end - beg);
        }

        for (uint16_t i_p = 1; i_p < p; i_p++) {
            blk_fps[i_p] = concat(blk_fps[i_p - 1], blk_fps[i_p], s * blks_per_thr);
        }

        fps.resize(num_blks + 1);
        no_init_resize(fps, num_blks);
        fps[0] = 0;

        #pragma omp parallel num_threads(p)
        {
            uint16_t i_p = omp_get_thread_num();
            const pos_t blk_beg = i_p * blks_per_thr;
            const pos_t blk_end = i_p == p - 1 ? num_blks : (blk_beg + blks_per_thr);
            fp_t fp = blk_fps[i_p];

            for (pos_t blk = blk_beg; blk < blk_end;) {
                pos_t beg = blk * s;
                pos_t end = beg + s;
                if (blk + 1 == num_blks) [[unlikely]] end = n;
                pos_t len = end - beg;
                fp = concat(fp, substring_naive<RIGHT>(beg, len), len);
                fps[++blk] = fp;
            }
        }
    }

    inline pos_t sampling_rate() const
    {
        return s;
    }

    inline pos_t window_size() const
    {
        return w;
    }

    inline uint64_t size_in_bytes() const
    {
        uint64_t size = sizeof(this);
        size += b_pow_leq_sqrt_n.size() * sizeof(fp_t);
        size += b_pow_step_sqrt_n.size() * sizeof(fp_t);
        size += fps.size() * sizeof(fp_t);
        return size;
    }

    template <typename char_t>
    inline fp_t push(const fp_t fp, const char_t chr) const
    {
        return mod(((b * fp_conc_t(fp)) + sqr_mers_prim) + fp_conc_t(char_to_uchar(chr)));
    }

    template <typename char_t>
    inline fp_t roll(const fp_t fp, const char_t pop, const char_t chr) const
    {
        return mod(((b * fp_conc_t(fp)) + pop_prec[char_to_uchar(pop)]) + fp_conc_t(char_to_uchar(chr)));
    }

    inline fp_t concat(const fp_t fp_lft, const fp_t fp_rght, const pos_t len_rght) const
    {
        return mod(fp_conc_t(b_pow(len_rght)) * fp_conc_t(fp_lft) + fp_conc_t(fp_rght));
    }

    template <direction dir = RIGHT>
    inline fp_t substring_naive(pos_t pos, const pos_t len) const
    {
        if constexpr (dir == LEFT) pos -= len - 1;
        fp_t fp = 0;

        for (pos_t i = pos; i < pos + len; i++) {
            fp = push(fp, T[i]);
        }

        return fp;
    }

    inline fp_t substring_up_to(const pos_t pos) const
    {
        pos_t blk = pos / s;
        pos_t blk_beg = blk * s;
        pos_t offs = pos - blk_beg;
        if (offs == 0) [[unlikely]] return fps[blk];
        return concat(fps[blk], substring_naive<RIGHT>(blk_beg, offs), offs);
    }

    template <direction dir = RIGHT>
    inline fp_t substring(pos_t pos, const pos_t len) const
    {
        if constexpr (dir == LEFT) pos -= len - 1;
        fp_t fp_beg_shft = mod(fp_conc_t(substring_up_to(pos)) * fp_conc_t(b_pow(len)));
        fp_t fp_end = substring_up_to(pos + len);
        return fp_end >= fp_beg_shft ? fp_end - fp_beg_shft : mers_prim - (fp_beg_shft - fp_end);
    }
};
