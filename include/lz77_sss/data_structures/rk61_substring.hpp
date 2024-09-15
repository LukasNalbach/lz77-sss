#pragma once

#include <cstdint>
#include <vector>
#include <random>

#include <lz77_sss/misc/utils.hpp>

class rk61_substring {
    protected:

    static constexpr uint64_t e61 = 61;
    static constexpr uint64_t p61 = (uint64_t(1) << e61) - 1;
    static constexpr uint128_t sq61 = uint128_t(p61) * p61;

    const char* T;
    uint64_t n;
    uint64_t sqrt_n;
    uint128_t b;
    uint64_t s;
    uint64_t w;
    uint128_t pop_prec[256] = {};
    std::vector<std::vector<uint64_t>> fp_smpl_tree;
    std::vector<uint64_t> bp_leq_sqrt_n;
    std::vector<uint64_t> bp_stp_sqrt_n;
    
    inline static uint64_t mod(const uint128_t val) {
        const uint128_t v = val + 1;
        const uint64_t z = ((v >> e61) + v) >> e61;
        return (val + z) & p61;
    }

    inline uint64_t b_pow(const uint64_t exp) const {
        const uint64_t exp_sqrt_n = exp / sqrt_n;
        const uint64_t offs_exp = exp - exp_sqrt_n * sqrt_n;
        return mod(uint128_t{bp_stp_sqrt_n[exp_sqrt_n]} *
                   uint128_t{bp_leq_sqrt_n[offs_exp]});
    }

    inline static uint64_t pow(uint64_t b, uint64_t exp) {
        uint64_t res = 1;
        while(exp > 0) {
            if(exp & 1ULL) {
                res = mod(uint128_t(b) * res);
            }
            b = mod(uint128_t(b) * b);
            exp >>= 1;
        }
        return res;
    }

    public:

    rk61_substring() = default;

    rk61_substring(
        const std::string& T,
        const uint64_t s,
        const uint64_t w = 0,
        uint16_t p = 1
    ) : T(T.data()), n(T.size()), s(s), w(w) {
        sqrt_n = std::ceil(std::sqrt(double(n)));
        std::random_device rd;
        std::mt19937_64 mt(rd());
        std::uniform_int_distribution<uint128_t> distrib(257, p61);
        b = distrib(mt);
        const uint64_t n_s = div_ceil<uint64_t>(n, s);

        no_init_resize(bp_leq_sqrt_n, sqrt_n + 1);
        no_init_resize(bp_stp_sqrt_n, sqrt_n + 1);
        bp_leq_sqrt_n[0] = 1;
        bp_stp_sqrt_n[0] = 1;

        for (uint64_t i = 1; i <= sqrt_n; i++) {
            bp_leq_sqrt_n[i] = mod(uint128_t{bp_leq_sqrt_n[i - 1]} * b);

            #ifndef NDEBUG
            assert(bp_leq_sqrt_n[i] == pow(b, i));
            #endif
        }

        for (uint64_t i = 1; i <= sqrt_n; i++) {
            bp_stp_sqrt_n[i] = mod(
                uint128_t{bp_stp_sqrt_n[i - 1]} *
                uint128_t{bp_leq_sqrt_n[sqrt_n]});
            
            #ifndef NDEBUG
            assert(bp_stp_sqrt_n[i] == pow(b, i * sqrt_n));
            #endif
        }

        if (w != 0) {
            const uint64_t max_exp_excl = b_pow(w);

            for(uint64_t c = 0; c < 256; c++) {
                pop_prec[c] = sq61 - uint128_t(max_exp_excl) * c;
            }
        }

        if (s >= n) return;
        const uint64_t h_max = log2_clz_up(n_s);
        fp_smpl_tree.resize(h_max + 1);
        no_init_resize(fp_smpl_tree[0], n_s);

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n_s; i++) {
            const uint64_t pos = i * s;
            const uint64_t len = std::min<uint64_t>((i + 1) * s, n) - pos;
            fp_smpl_tree[0][i] = substring_naive(pos, len);
        }

        for (uint64_t h = 1; h <= h_max; h++) {
            const uint64_t ch = h - 1;
            const uint64_t w_h = div_ceil<uint64_t>(n_s, 1 << h);
            const uint64_t l_h = s * (1 << h);
            const uint64_t l_ch = s * (1 << ch);
            no_init_resize(fp_smpl_tree[h], w_h);

            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < w_h - 1; i++) {
                fp_smpl_tree[h][i] = concat(
                    fp_smpl_tree[ch][2 * i],
                    fp_smpl_tree[ch][2 * i + 1],
                    l_ch
                );

                #ifndef NDEBUG
                assert(fp_smpl_tree[h][i] == substring_naive(i * l_h, l_h));
                #endif
            }

            uint64_t p_lc = 2 * (w_h - 1);
            uint64_t p_rc = 2 * (w_h - 1) + 1;

            if (p_rc * l_ch >= n) {
                fp_smpl_tree[h][w_h - 1] = fp_smpl_tree[ch][p_lc];
            } else {
                fp_smpl_tree[h][w_h - 1] = concat(
                    fp_smpl_tree[ch][p_lc],
                    fp_smpl_tree[ch][p_rc],
                    n - p_rc * l_ch
                );

                #ifndef NDEBUG
                assert(fp_smpl_tree[h][w_h - 1] ==
                    substring_naive((w_h - 1) * l_h, n - p_lc * l_ch));
                #endif
            }
        }
    }

    inline uint64_t sampling_rate() const {
        return s;
    }

    inline uint64_t window_size() const {
        return w;
    }

    inline uint64_t size_in_bytes() const {
        uint64_t size = sizeof(this);
        size += bp_leq_sqrt_n.size() * sizeof(uint64_t);
        size += bp_stp_sqrt_n.size() * sizeof(uint64_t);
        
        for (auto& vec : fp_smpl_tree) {
            size += vec.size() * sizeof(uint64_t);
        }

        return size;
    }

    inline uint64_t push(const uint64_t fp, const char push) const {
        return mod(((b * uint128_t{fp}) + sq61) + uint128_t(char_to_uchar(push)));
    }

    inline uint64_t roll(const uint64_t fp, const char pop, const char push) const {
        return mod(((b * uint128_t{fp}) + pop_prec[char_to_uchar(pop)]) + uint128_t(char_to_uchar(push)));
    }

    inline uint64_t concat(const uint64_t fp_l, const uint64_t fp_r, const uint64_t len_r) const {
        return mod(uint128_t(b_pow(len_r)) * uint128_t(fp_l) + uint128_t(fp_r));
    }

    template <direction dir = RIGHT>
    inline uint64_t substring_naive(uint64_t pos, const uint64_t len) const {
        if constexpr (dir == LEFT) pos -= len - 1;
        uint64_t fp = 0;

        for (uint64_t i = pos; i < pos + len; i++) {
            fp = push(fp, T[i]);
        }

        return fp;
    }

    template <direction dir = RIGHT>
    inline uint64_t substring(uint64_t pos, const uint64_t len) const {
        if constexpr (dir == LEFT) {
            pos -= len - 1;
        }

        if (len < s) [[likely]] {
            return substring_naive<RIGHT>(pos, len);
        }

        uint64_t h = log2_clz(len / s) - 1;
        uint64_t sl_h = s * (1 << h);
        uint64_t i = div_ceil<uint64_t>(pos, sl_h);
        const uint64_t end = pos + len;

        if ((i + 1) * sl_h > end) {
            if (h == 0) return substring_naive<RIGHT>(pos, len);
            h--;
            sl_h /= 2;
            i = div_ceil<uint64_t>(pos, sl_h);
        }

        const uint64_t p_m = i * sl_h;
        const uint64_t p_r = p_m + sl_h;
        uint64_t fp = fp_smpl_tree[h][i];
        
        if (p_m > pos) [[likely]] {
            const uint64_t fp_l = substring<RIGHT>(pos, p_m - pos);
            fp = concat(fp_l, fp, sl_h);
        }
        
        if (p_r < end) [[likely]] {
            const uint64_t l_r = end - p_r;
            const uint64_t fp_r = substring<RIGHT>(p_r, l_r);
            fp = concat(fp, fp_r, l_r);
        }
        
        return fp;
    }
};
