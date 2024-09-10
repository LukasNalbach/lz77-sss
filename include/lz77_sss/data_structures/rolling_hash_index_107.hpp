#pragma once

#include <lz77_sss/misc/utils.hpp>

template <typename pos_t, uint8_t num_patt_lens>
class rolling_hash_index_107 {
    protected:

    using rolling_hash_t = alx::rolling_hash::rk_prime<107>;
    using fingerprint_t = uint128_t;
    char* input = nullptr;
    pos_t input_size = 0;
    std::array<pos_t, num_patt_lens> patt_lens;
    rolling_hash_t* rolling_hash[num_patt_lens] = {nullptr};
    std::vector<pos_t> H;
    fingerprint_t h_mod_mask = 0;
    pos_t cur_pos = 0;

    public:
    rolling_hash_index_107() = default;

    rolling_hash_index_107(
        char* input, pos_t size,
        std::array<pos_t, num_patt_lens> patt_lens,
        uint64_t target_size_in_bytes
    ) : input(input), input_size(size), patt_lens(patt_lens) {
        uint64_t target_size_h = std::max<int64_t>(
            input_size / (sizeof(pos_t) * 10),
            (int64_t{target_size_in_bytes} - int64_t{sizeof(rolling_hash_t) *
            num_patt_lens}) /int64_t{sizeof(pos_t)});
        uint8_t log2_size_h = std::round(std::log2(target_size_h));
        pos_t size_h = 1 << log2_size_h;
        h_mod_mask = size_h - 1;
        no_init_resize(H, size_h);
        std::fill(H.begin(), H.end(), std::numeric_limits<pos_t>::max());

        for_constexpr<0, num_patt_lens, 1>([&](auto i) {
            rolling_hash[i] = new rolling_hash_t(
                fingerprint_t{patt_lens[i]});
        });

        reinit(0);
    }

    inline void reinit(pos_t pos) {
        cur_pos = pos;

        for_constexpr<0, num_patt_lens, 1>([&](auto i) {
            rolling_hash[i]->reset();
            init<i>();
        });
    }
    
    template <pos_t i>
    inline void init() {
        for (pos_t j = 0; j < patt_lens[i]; j++) {
            rolling_hash[i]->roll_in(
                char_to_uchar(input[cur_pos + j]));
        }
    }

    template <pos_t i>
    inline void roll() {
        rolling_hash[i]->roll(
            char_to_uchar(input[cur_pos]),
            char_to_uchar(input[cur_pos + patt_lens[i]])
        );
    }

    inline void roll() {
        for_constexpr<0, num_patt_lens, 1>([&](auto i) {
            roll<i>();
        });

        cur_pos++;
    }

    template <pos_t i>
    inline void advance() {
        if (cur_pos + patt_lens[i] < input_size) {
            H[rolling_hash[i]->get_fp() & h_mod_mask] = cur_pos;
            roll<i>();
        }
    }

    inline void advance() {
        for_constexpr<0, num_patt_lens, 1>([&](auto i) {
            advance<i>();
        });

        cur_pos++;
    }

    template <pos_t i>
    inline pos_t advance_and_get_occ() {
        fingerprint_t h_pos = rolling_hash[i]->get_fp() & h_mod_mask;
        pos_t occ = H[h_pos];
        H[h_pos] = cur_pos;

        if (cur_pos + patt_lens[i] < input_size) {
            roll<i>();
        }

        return occ;
    }

    template <pos_t i>
    inline pos_t get_occ() const {
        return H[rolling_hash[i]->get_fp() & h_mod_mask];
    }

    inline uint64_t size_in_bytes() const {
        return sizeof(pos_t) * H.size() + sizeof(this) +
            sizeof(rolling_hash_t) * num_patt_lens;
    }

    double rate_init() const {
        pos_t num_init = 0;

        for (pos_t i = 0; i < H.size(); i++) {
            if (H[i] < input_size) {
                num_init++;
            }
        }

        return num_init / (double) H.size();
    }

    inline pos_t pos() const {
        return cur_pos;
    }

    inline void inc_pos() {
        cur_pos++;
    }
};