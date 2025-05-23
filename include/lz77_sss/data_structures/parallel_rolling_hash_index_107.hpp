#pragma once

template <typename pos_t> class lz77_sss;

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t, uint8_t num_patt_lens, typename char_t>
class parallel_rolling_hash_index_107 {
public:
    using fp_t = uint128_t;
    using fp_arr_t = std::array<fp_t, num_patt_lens>;

protected:
    using rolling_hash_t = lce::rolling_hash::rk_prime<107>;
    const char_t* input = nullptr;
    pos_t input_size = 0;
    std::array<pos_t, num_patt_lens> patt_lens;
    rolling_hash_t* rolling_hash[num_patt_lens] = { nullptr };
    std::vector<pos_t> H_old;
    std::vector<pos_t> H_new;
    fp_t h_mod_mask = 0;

public:
    parallel_rolling_hash_index_107() = default;

    parallel_rolling_hash_index_107(
        const char_t* input, pos_t size,
        std::array<pos_t, num_patt_lens> patt_lens,
        int64_t target_size_in_bytes, uint16_t num_threads)
        : input(input)
        , input_size(size)
        , patt_lens(patt_lens)
    {
        static constexpr int64_t rolling_hash_size = sizeof(rolling_hash_t) * num_patt_lens;
        int64_t min_index_size = std::max<pos_t>(lz77_sss<pos_t>::min_rh_index_size,
            input_size * lz77_sss<pos_t>::min_rel_rh_index_size) / sizeof(pos_t);
        int64_t max_index_size = lz77_sss<pos_t>::max_rh_index_size / sizeof(pos_t);
        int64_t target_index_size = std::max<int64_t>(0, target_size_in_bytes - rolling_hash_size) / sizeof(pos_t);

        uint64_t target_size_h = std::min<int64_t>(max_index_size, std::max<int64_t>(min_index_size, target_index_size));
        uint8_t log2_size_h = std::round(std::log2(target_size_h)) - 1;
        pos_t size_h = 1 << log2_size_h;
        h_mod_mask = size_h - 1;
        no_init_resize(H_old, size_h);
        no_init_resize(H_new, size_h);

        uint64_t num_blks = H_new.size() / (sizeof(uint128_t) / sizeof(pos_t));
        uint128_t* blks_new = reinterpret_cast<uint128_t*>(H_new.data());

        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t i = 0; i < num_blks; i++) {
            blks_new[i] = std::numeric_limits<uint128_t>::max();
        }

        for_constexpr<0, num_patt_lens, 1>([&](auto i) {
            rolling_hash[i] = new rolling_hash_t(
                fp_t { patt_lens[i] });
        });
    }

    inline uint64_t size() const
    {
        return H_old.size();
    }

    void overwrite(uint16_t p)
    {
        uint128_t* blks_old = reinterpret_cast<uint128_t*>(H_old.data());
        uint128_t* blks_new = reinterpret_cast<uint128_t*>(H_new.data());
        uint64_t num_blks = H_new.size() / (sizeof(uint128_t) / sizeof(pos_t));

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < num_blks; i++) {
            blks_old[i] = blks_new[i];
        }
    }

    inline void reinit(fp_arr_t& fps, pos_t pos)
    {
        for_constexpr<0, num_patt_lens, 1>([&](auto i) {
            fps[i] = 0;

            for (pos_t j = 0; j < patt_lens[i]; j++) {
                if (pos + patt_lens[i] < input_size) [[likely]] {
                    fps[i] = rolling_hash[i]->roll_in(
                        fps[i], char_to_uchar(input[pos + j]));
                }
            }
        });
    }

    template <pos_t i>
    inline void roll(fp_arr_t& fps, pos_t pos)
    {
        fps[i] = rolling_hash[i]->roll(
            fps[i],
            char_to_uchar(input[pos]),
            char_to_uchar(input[pos + patt_lens[i]]));
    }

    inline void roll(fp_arr_t& fps, pos_t pos)
    {
        for_constexpr<0, num_patt_lens, 1>([&](auto i) {
            if (pos + patt_lens[i] < input_size) [[likely]] {
                roll<i>(fps, pos);
            }
        });
    }

    template <pos_t i>
    inline void advance(fp_arr_t& fps, pos_t pos)
    {
        if (pos + patt_lens[i] < input_size) [[likely]] {
            H_new[fps[i] & h_mod_mask] = pos;
            roll<i>(fps, pos);
        }
    }

    inline void advance(fp_arr_t& fps, pos_t pos)
    {
        for_constexpr<0, num_patt_lens, 1>([&](auto i) {
            advance<i>(fps, pos);
        });
    }

    template <pos_t i, bool ret_new = false>
    inline pos_t advance_and_get_occ(fp_arr_t& fps, pos_t pos)
    {
        fp_t h_pos = fps[i] & h_mod_mask;
        pos_t occ;

        if constexpr (ret_new) {
            occ = H_new[h_pos];
        }
        
        H_new[h_pos] = pos;

        if (pos + patt_lens[i] < input_size) [[likely]] {
            roll<i>(fps, pos);
        }

        if constexpr (ret_new) {
            return occ;
        } else {
            return H_old[h_pos];
        }
    }

    template <pos_t i>
    inline pos_t get_occ(fp_arr_t& fps) const
    {
        return H_old[fps[i] & h_mod_mask];
    }

    inline uint64_t size_in_bytes() const
    {
        return sizeof(this) + 2 * sizeof(pos_t) * H_old.size() +
            sizeof(rolling_hash_t) * num_patt_lens;
    }
};