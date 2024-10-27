#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t, uint8_t num_patt_lens>
class parallel_rolling_hash_index_107 {
public:
    using fp_t = uint128_t;
    using fp_arr_t = std::array<fp_t, num_patt_lens>;

protected:
    using rolling_hash_t = alx::rolling_hash::rk_prime<107>;
    static constexpr double min_rel_idx_size = 0.1;
    char* input = nullptr;
    pos_t input_size = 0;
    std::array<pos_t, num_patt_lens> patt_lens;
    rolling_hash_t* rolling_hash[num_patt_lens] = { nullptr };
    std::vector<pos_t> H_old;
    std::vector<pos_t> H_new;
    fp_t h_mod_mask = 0;

public:
    parallel_rolling_hash_index_107() = default;

    parallel_rolling_hash_index_107(
        char* input, pos_t size,
        std::array<pos_t, num_patt_lens> patt_lens,
        uint64_t target_size_in_bytes, uint16_t num_threads)
        : input(input)
        , input_size(size)
        , patt_lens(patt_lens)
    {
        uint64_t target_size_h = std::max<int64_t>(
            (input_size * min_rel_idx_size) / (2 * sizeof(pos_t)),
            (int64_t { target_size_in_bytes } - int64_t { sizeof(rolling_hash_t) *
            num_patt_lens }) / int64_t { 2 * sizeof(pos_t) });
        uint8_t log2_size_h = std::round(std::log2(target_size_h));
        pos_t size_h = 1 << log2_size_h;
        h_mod_mask = size_h - 1;
        H_old.resize(size_h, { });
        H_new.resize(size_h, { });
        uint64_t num_elements = H_old.size() / (16 / sizeof(pos_t));
        uint128_t* data_new = reinterpret_cast<uint128_t*>(H_new.data());

        for_constexpr<0, num_patt_lens, 1>([&](auto i) {
            rolling_hash[i] = new rolling_hash_t(
                fp_t { patt_lens[i] });
        });

        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t i = 0; i < num_elements; i++) {
            data_new[i] = std::numeric_limits<uint128_t>::max();
        }
    }

    inline uint64_t size() const
    {
        return H_old.size();
    }

    inline uint128_t* data_old()
    {
        return reinterpret_cast<uint128_t*>(H_old.data());
    }

    inline uint128_t* data_new()
    {
        return reinterpret_cast<uint128_t*>(H_new.data());
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
            roll<i>(fps, pos);
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
        H_new[h_pos] = pos;

        if (pos + patt_lens[i] < input_size) [[likely]] {
            roll<i>(fps, pos);
        }

        if constexpr (ret_new) {
            return H_new[h_pos];
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