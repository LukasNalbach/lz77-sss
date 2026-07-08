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

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <istream>
#include <ostream>
#include <queue>
#include <span>
#include <vector>

__extension__ typedef unsigned __int128 uint128_t;

class bit_writer {
    std::ostream& out;
    std::vector<char> chunk;
    uint64_t chunk_pos = 0;
    uint128_t acc = 0;
    uint64_t cnt = 0;
    static constexpr uint64_t chunk_size = 1 << 16;

    inline void put_byte(char b) {
        chunk[chunk_pos++] = b;
        if (chunk_pos == chunk_size) {
            out.write(chunk.data(), chunk_size);
            chunk_pos = 0;
        }
    }

public:
    explicit bit_writer(std::ostream& out) : out(out), chunk(chunk_size) {}

    inline void put(uint64_t x, uint64_t n) {
        if (n == 0) return;
        
        if (chunk_pos > chunk_size - 8) {
            out.write(chunk.data(), chunk_pos);
            chunk_pos = 0;
        }

        uint64_t mask = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
        acc |= uint128_t{x & mask} << (128 - cnt - n);
        cnt += n;

        if (cnt >= 64) {
            uint64_t to_write = uint64_t(acc >> 64);
            uint64_t sw = std::byteswap(to_write);
            std::memcpy(&chunk[chunk_pos], &sw, 8);
            chunk_pos += 8;
            acc <<= 64;
            cnt -= 64;
        }
    }

    inline void flush() {
        uint64_t bytes = (cnt + 7) / 8;
        for (uint64_t i = 0; i < bytes; ++i) {
            put_byte((char)(acc >> 120));
            acc <<= 8;
        }
        cnt = 0;
        acc = 0;
        if (chunk_pos > 0) {
            out.write(chunk.data(), chunk_pos);
            chunk_pos = 0;
        }
    }
};

class bit_reader {
    std::istream& in;
    std::vector<char> chunk;
    uint64_t chunk_pos = 0;
    uint64_t chunk_len = 0;
    uint128_t acc = 0;
    uint64_t cnt = 0;
    static constexpr uint64_t chunk_size = 1 << 16;

    inline void refill() {
        if (chunk_pos + 8 <= chunk_len) {
            uint64_t raw;
            std::memcpy(&raw, &chunk[chunk_pos], 8);
            uint64_t next_val = std::byteswap(raw);
            acc |= uint128_t{next_val} << (64 - cnt);
            cnt += 64;
            chunk_pos += 8;
        } else {
            while (cnt <= 120 && chunk_pos < chunk_len) {
                uint64_t b = (uint8_t)(chunk[chunk_pos++]);
                acc |= uint128_t{b} << (128 - cnt - 8);
                cnt += 8;
            }
        }
    }

public:
    explicit bit_reader(std::istream& in) : in(in), chunk(chunk_size) {}

    template <typename T>
    inline T get_bit() {
        if (cnt < 1) {
            if (chunk_pos >= chunk_len) {
                in.read(chunk.data(), chunk_size);
                chunk_len = in.gcount();
                chunk_pos = 0;
            }
            refill();
        }
        T bit = acc >> 127;
        acc <<= 1;
        cnt--;
        return bit;
    }

    inline uint64_t get(uint64_t n) {
        if (n == 0) return 0;
        if (cnt < n) {
            if (chunk_pos >= chunk_len) {
                in.read(chunk.data(), chunk_size);
                chunk_len = in.gcount();
                chunk_pos = 0;
            }
            refill();
        }
        uint64_t res = acc >> (128 - n);
        acc <<= n;
        cnt -= n;
        return res;
    }
};

inline void put_elias_delta(bit_writer& out, uint64_t x) {
    if (x == 0) return;
    uint64_t len_x = std::bit_width(x);
    uint64_t len_len = std::bit_width(len_x);
    out.put(0, len_len - 1);
    out.put(len_x, len_len);
    out.put(x, len_x - 1);
}

inline uint64_t get_elias_delta(bit_reader& in) {
    uint64_t zeros = 0;
    while (in.get_bit<uint8_t>() == 0) zeros++;
    uint64_t len_x = (1ULL << zeros) | in.get(zeros);
    if (len_x <= 1) return len_x;
    return (1ULL << (len_x - 1)) | in.get(len_x - 1);
}

class huffman {
public:
    static constexpr uint64_t max_len = 15;

private:
    uint64_t sigma = 0;
    std::vector<uint64_t> length; 
    std::vector<uint64_t> code;
    std::array<uint64_t, max_len + 1> cnt{};
    std::array<uint64_t, max_len + 1> first_code{};
    std::array<uint64_t, max_len + 1> first_index{};
    std::vector<uint64_t> sorted_syms;

    void build_codes() {
        cnt.fill(0);
        for (uint64_t s = 0; s < sigma; s++) {
            if (length[s]) cnt[length[s]]++;
        }

        code.assign(sigma, 0);
        std::array<uint64_t, max_len + 1> next{};
        uint64_t c = 0;
        for (uint64_t l = 1; l <= max_len; l++) { 
            c = (c + cnt[l - 1]) << 1; 
            next[l] = c; 
        }
        for (uint64_t s = 0; s < sigma; s++) {
            if (length[s]) code[s] = next[length[s]]++;
        }

        sorted_syms.clear();
        for (uint64_t l = 1; l <= max_len; l++) {
            for (uint64_t s = 0; s < sigma; s++) {
                if (length[s] == l) sorted_syms.push_back(s);
            }
        }

        uint64_t cc = 0;
        uint64_t idx = 0;
        for (uint64_t l = 1; l <= max_len; l++) {
            first_code[l] = cc; 
            first_index[l] = idx;
            cc = (cc + cnt[l]) << 1; 
            idx += cnt[l];
        }
    }

public:
    void build_from_freq(std::span<const uint64_t> freq) {
        sigma = freq.size();
        length.assign(sigma, 0);

        std::vector<uint64_t> used;
        for (uint64_t s = 0; s < sigma; s++) {
            if (freq[s]) used.push_back(s);
        }

        if (used.empty()) { build_codes(); return; }
        if (used.size() == 1) { length[used[0]] = 1; build_codes(); return; }

        static constexpr uint64_t no_idx = ~0ULL;
        struct node { uint64_t w; uint64_t l, r, sym; };
        std::vector<node> nodes;
        nodes.reserve(2 * used.size());
        
        using item = std::pair<uint64_t, uint64_t>;
        std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
        
        for (uint64_t s : used) {
            nodes.push_back({ freq[s], no_idx, no_idx, s });
            pq.push({ freq[s], nodes.size() - 1 });
        }
        
        while (pq.size() > 1) {
            auto a = pq.top(); pq.pop();
            auto b = pq.top(); pq.pop();
            nodes.push_back({ a.first + b.first, a.second, b.second, no_idx });
            pq.push({ a.first + b.first, nodes.size() - 1 });
        }

        std::vector<uint64_t> nat(sigma, 0);
        std::vector<std::pair<uint64_t, uint64_t>> st { { pq.top().second, 0 } };
        while (!st.empty()) {
            auto [u, d] = st.back(); st.pop_back();
            if (nodes[u].sym != no_idx) nat[nodes[u].sym] = std::max<uint64_t>(1, d);
            else { st.push_back({ nodes[u].l, d + 1 }); st.push_back({ nodes[u].r, d + 1 }); }
        }

        std::array<uint64_t, max_len + 2> bl{};
        for (uint64_t s : used) {
            bl[std::min<uint64_t>(nat[s], max_len)]++;
        }
        
        uint64_t kraft = 0;
        for (uint64_t l = 1; l <= max_len; l++) kraft += bl[l] << (max_len - l);
        
        uint64_t target = (1ULL << max_len);
        while (kraft > target) {
            uint64_t l = max_len - 1;
            while (l >= 1 && bl[l] == 0) l--;
            bl[l]--; bl[l + 1]++;
            kraft -= (1ULL << (max_len - l - 1));
        }

        std::sort(used.begin(), used.end(), [&](uint64_t a, uint64_t b) { return freq[a] < freq[b]; });
        uint64_t idx = 0;
        for (uint64_t l = max_len; l >= 1; l--) {
            for (uint64_t k = 0; k < bl[l]; k++) {
                length[used[idx++]] = l;
            }
        }

        build_codes();
    }

    void write_table(bit_writer& out) const {
        for (uint64_t s = 0; s < sigma; s++) out.put(length[s], 4);
    }

    void read_table(bit_reader& in, uint64_t alphabet) {
        sigma = alphabet;
        length.assign(sigma, 0);
        for (uint64_t s = 0; s < sigma; s++) length[s] = in.get(4);
        build_codes();
    }

    inline void encode(bit_writer& out, uint64_t sym) const {
        out.put(code[sym], length[sym]);
    }

    inline uint64_t decode(bit_reader& in) const {
        uint64_t c = 0;
        for (uint64_t l = 1; l <= max_len; l++) {
            c = (c << 1) | in.get_bit<uint64_t>();
            if (cnt[l] && c - first_code[l] < cnt[l])
                return sorted_syms[first_index[l] + (c - first_code[l])];
        }
        return 0;
    }
};

static constexpr uint32_t huff_sigma = 66;
static constexpr uint64_t huff_block_size = 1 << 14;

class huff_writer {
    bit_writer out;
    struct entry { uint64_t val; uint64_t len; };
    std::vector<entry> block;
    uint64_t pos = 0;

    void flush_block() {
        if (block.empty()) return;

        std::array<uint64_t, huff_sigma> hl{}, hd{};
        for (const entry& e : block) {
            hl[e.len == 0 ? 0 : std::bit_width(e.len)]++;
            if (e.len != 0) hd[std::bit_width(e.val)]++;
        }

        huffman len_huff, dist_huff;
        len_huff.build_from_freq(hl);
        dist_huff.build_from_freq(hd);

        put_elias_delta(out, block.size());
        len_huff.write_table(out);
        dist_huff.write_table(out);

        for (const entry& e : block) {
            uint64_t lb = e.len == 0 ? 0 : std::bit_width(e.len);
            len_huff.encode(out, lb);
            if (lb == 0) {
                out.put(e.val, 8);
            } else {
                out.put(e.len, lb - 1);
                uint64_t db = std::bit_width(e.val);
                dist_huff.encode(out, db);
                out.put(e.val, db - 1);
            }
        }

        block.clear();
    }

public:
    huff_writer(std::ostream& o, uint64_t n) : out(o) {
        for (int i = 0; i < 5; ++i) {
            o.put(n);
            n >>= 8;
        }
        block.reserve(huff_block_size);
        pos = 0;
    }

    template <typename factor_t>
    void add(const factor_t& f) {
        if (f.len == 0) {
            block.push_back({ f.src & 0xFFu, 0 });
            pos += 1;
        } else {
            block.push_back({ pos - f.src, f.len });
            pos += f.len;
        }
        if (block.size() == huff_block_size) flush_block();
    }

    void finish() {
        flush_block();
        out.flush();
    }
};

template <typename factor_t>
class huff_factor_iterator {
    bit_reader& r;
    huffman& len_huff;
    huffman& dist_huff;
    uint64_t remaining = 0;
    uint64_t pos = 0;
    uint64_t n = 0;
    factor_t cur{};

    void advance() {
        if (pos >= n) return;

        if (remaining == 0) {
            remaining = get_elias_delta(r);
            len_huff.read_table(r, huff_sigma);
            dist_huff.read_table(r, huff_sigma);
        }

        uint64_t lb = len_huff.decode(r);
        if (lb == 0) {
            cur.src = r.get(8);
            cur.len = 0;
            pos += 1;
        } else {
            uint64_t len = (1ULL << (lb - 1)) | r.get(lb - 1);
            uint64_t db = dist_huff.decode(r);
            uint64_t dist = (1ULL << (db - 1)) | r.get(db - 1);
            cur.len = len;
            cur.src = pos - dist;
            pos += len;
        }
        remaining--;
    }

public:
    using iterator_category = std::input_iterator_tag;
    using value_type = factor_t;
    using difference_type = std::ptrdiff_t;
    using pointer = const factor_t*;
    using reference = const factor_t&;

    huff_factor_iterator() = default;
    huff_factor_iterator(bit_reader& reader, huffman& lh, huffman& dh, uint64_t n)
        : r(reader), len_huff(lh), dist_huff(dh), n(n) { advance(); }

    reference operator*() const { return cur; }
    huff_factor_iterator& operator++() { advance(); return *this; }
    huff_factor_iterator operator++(int) { auto tmp = *this; advance(); return tmp; }
};