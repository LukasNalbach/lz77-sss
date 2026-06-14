#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <istream>
#include <ostream>
#include <queue>
#include <string>
#include <vector>

class bit_writer {
    std::ostream& out;
    std::string chunk;
    size_t chunk_pos = 0;
    uint64_t acc = 0;
    uint32_t cnt = 0;
    static constexpr size_t chunk_size = 1 << 16;

    inline void put_byte(uint8_t b)
    {
        chunk[chunk_pos++] = (char) b;
        if (chunk_pos == chunk_size) {
            out.write(chunk.data(), chunk_size);
            chunk_pos = 0;
        }
    }

public:
    bit_writer(std::ostream& out) : out(out), chunk(chunk_size, '\0') { }

    inline void put(uint64_t x, uint32_t n)
    {
        if (n > 32) {
            put(x >> 32, n - 32);
            n = 32;
        }

        uint32_t low = (n == 32) ? (uint32_t) x : (uint32_t) (x & (((uint32_t) 1 << n) - 1));
        acc = (acc << n) | low;
        cnt += n;

        while (cnt >= 8) {
            cnt -= 8;
            put_byte((uint8_t) (acc >> cnt));
        }
    }

    inline void flush()
    {
        if (cnt > 0) {
            put_byte((uint8_t) (acc << (8 - cnt)));
            acc = 0;
            cnt = 0;
        }
        if (chunk_pos > 0) {
            out.write(chunk.data(), chunk_pos);
            chunk_pos = 0;
        }
    }
};

class bit_reader {
    std::istream& in;
    std::string chunk;
    size_t chunk_pos = 0;
    size_t chunk_len = 0;
    uint64_t acc = 0;
    uint32_t cnt = 0;
    static constexpr size_t chunk_size = 1 << 16;

    inline uint8_t next_byte()
    {
        if (chunk_pos >= chunk_len) {
            in.read(chunk.data(), chunk_size);
            chunk_len = (size_t) in.gcount();
            chunk_pos = 0;
            if (chunk_len == 0) return 0;
        }
        return (uint8_t) chunk[chunk_pos++];
    }

public:
    bit_reader(std::istream& in) : in(in), chunk(chunk_size, '\0') { }

    inline uint32_t get_bit()
    {
        if (cnt == 0) {
            acc = (acc << 8) | next_byte();
            cnt = 8;
        }
        cnt--;
        return (uint32_t) ((acc >> cnt) & 1);
    }

    inline uint64_t get(uint32_t n)
    {
        if (n > 32) {
            uint64_t hi = get(n - 32);
            return (hi << 32) | get(32);
        }

        while (cnt < n) {
            acc = (acc << 8) | next_byte();
            cnt += 8;
        }

        cnt -= n;
        uint32_t mask = (n == 32) ? ~0u : (((uint32_t) 1 << n) - 1);
        return (acc >> cnt) & mask;
    }
};

inline void put_elias_delta(bit_writer& out, uint64_t x)
{
    uint32_t len_x = std::bit_width(x);
    uint32_t len_len = std::bit_width(len_x);
    out.put(0, len_len - 1);
    out.put(len_x, len_len);
    out.put(x, len_x - 1);
}

inline uint64_t get_elias_delta(bit_reader& in)
{
    uint32_t zeros = 0;
    while (in.get_bit() == 0) zeros++;
    uint32_t len_x = ((uint32_t) 1 << zeros) | (uint32_t) in.get(zeros);
    return ((uint64_t) 1 << (len_x - 1)) | in.get(len_x - 1);
}

class huffman {
public:
    static constexpr uint8_t max_len = 15;

private:
    uint32_t sigma = 0;
    std::vector<uint8_t> length;
    std::vector<uint32_t> code;
    std::array<uint32_t, max_len + 1> cnt {};
    std::array<uint32_t, max_len + 1> first_code {};
    std::array<uint32_t, max_len + 1> first_index {};
    std::vector<uint32_t> sorted_syms;

    void build_codes()
    {
        cnt.fill(0);
        for (uint32_t s = 0; s < sigma; s++)
            if (length[s]) cnt[length[s]]++;

        code.assign(sigma, 0);
        std::array<uint32_t, max_len + 1> next {};
        uint32_t c = 0;
        for (uint8_t l = 1; l <= max_len; l++) { c = (c + cnt[l - 1]) << 1; next[l] = c; }
        for (uint32_t s = 0; s < sigma; s++)
            if (length[s]) code[s] = next[length[s]]++;

        sorted_syms.clear();
        for (uint8_t l = 1; l <= max_len; l++)
            for (uint32_t s = 0; s < sigma; s++)
                if (length[s] == l) sorted_syms.push_back(s);

        uint32_t cc = 0, idx = 0;
        for (uint8_t l = 1; l <= max_len; l++) {
            first_code[l] = cc; first_index[l] = idx;
            cc = (cc + cnt[l]) << 1; idx += cnt[l];
        }
    }

public:
    void build_from_freq(const std::vector<uint64_t>& freq)
    {
        sigma = (uint32_t) freq.size();
        length.assign(sigma, 0);

        std::vector<uint32_t> used;
        for (uint32_t s = 0; s < sigma; s++) if (freq[s]) used.push_back(s);

        if (used.empty()) { build_codes(); return; }
        if (used.size() == 1) { length[used[0]] = 1; build_codes(); return; }

        struct node { uint64_t w; int l, r, sym; };
        std::vector<node> nodes;
        nodes.reserve(2 * used.size());
        using item = std::pair<uint64_t, int>;
        std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
        for (uint32_t s : used) {
            nodes.push_back({ freq[s], -1, -1, (int) s });
            pq.push({ freq[s], (int) nodes.size() - 1 });
        }
        while (pq.size() > 1) {
            auto a = pq.top(); pq.pop();
            auto b = pq.top(); pq.pop();
            nodes.push_back({ a.first + b.first, a.second, b.second, -1 });
            pq.push({ a.first + b.first, (int) nodes.size() - 1 });
        }

        std::vector<uint32_t> nat(sigma, 0);
        std::vector<std::pair<int, uint32_t>> st { { pq.top().second, 0 } };
        while (!st.empty()) {
            auto [u, d] = st.back(); st.pop_back();
            if (nodes[u].sym >= 0) nat[nodes[u].sym] = std::max<uint32_t>(1, d);
            else { st.push_back({ nodes[u].l, d + 1 }); st.push_back({ nodes[u].r, d + 1 }); }
        }

        std::vector<uint32_t> bl(max_len + 2, 0);
        for (uint32_t s : used) bl[std::min<uint32_t>(nat[s], max_len)]++;
        uint64_t kraft = 0;
        for (uint32_t l = 1; l <= max_len; l++) kraft += (uint64_t) bl[l] << (max_len - l);
        uint64_t target = (uint64_t) 1 << max_len;
        while (kraft > target) {
            int l = max_len - 1;
            while (l >= 1 && bl[l] == 0) l--;
            bl[l]--; bl[l + 1]++;
            kraft -= (uint64_t) 1 << (max_len - l - 1);
        }

        std::sort(used.begin(), used.end(), [&](uint32_t a, uint32_t b) { return freq[a] < freq[b]; });
        uint32_t idx = 0;
        for (int l = max_len; l >= 1; l--)
            for (uint32_t k = 0; k < bl[l]; k++)
                length[used[idx++]] = (uint8_t) l;

        build_codes();
    }

    void write_table(bit_writer& out) const
    {
        for (uint32_t s = 0; s < sigma; s++) out.put(length[s], 4);
    }

    void read_table(bit_reader& in, uint32_t alphabet)
    {
        sigma = alphabet;
        length.assign(sigma, 0);
        for (uint32_t s = 0; s < sigma; s++) length[s] = (uint8_t) in.get(4);
        build_codes();
    }

    inline void encode(bit_writer& out, uint32_t sym) const
    {
        out.put(code[sym], length[sym]);
    }

    inline uint32_t decode(bit_reader& in) const
    {
        uint32_t c = 0;
        for (uint8_t l = 1; l <= max_len; l++) {
            c = (c << 1) | in.get_bit();
            if (cnt[l] && c - first_code[l] < cnt[l])
                return sorted_syms[first_index[l] + (c - first_code[l])];
        }
        return 0;
    }
};

static constexpr uint32_t huff_sigma = 66;
static constexpr size_t huff_block_size = 1 << 14;

class huff_writer {
    bit_writer out;
    struct entry { uint64_t val; uint64_t len; };
    std::vector<entry> block;
    uint64_t pos = 0;

    void flush_block()
    {
        if (block.empty()) return;

        std::vector<uint64_t> hl(huff_sigma, 0), hd(huff_sigma, 0);
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
            uint8_t lb = e.len == 0 ? 0 : (uint8_t) std::bit_width(e.len);
            len_huff.encode(out, lb);
            if (lb == 0) {
                out.put(e.val, 8);
            } else {
                out.put(e.len, lb - 1);
                uint8_t db = (uint8_t) std::bit_width(e.val);
                dist_huff.encode(out, db);
                out.put(e.val, db - 1);
            }
        }

        block.clear();
    }

public:
    huff_writer(std::ostream& o, uint64_t n) : out(o)
    {
        o.write((char*) &n, 5);
        block.reserve(huff_block_size);
    }

    template <typename factor_t>
    void add(const factor_t& f)
    {
        if (f.len == 0) {
            block.push_back({ (uint64_t) (uint8_t) f.src, 0 });
            pos++;
        } else {
            block.push_back({ pos - (uint64_t) f.src, (uint64_t) f.len });
            pos += f.len;
        }
        if (block.size() == huff_block_size) flush_block();
    }

    void finish()
    {
        flush_block();
        out.flush();
    }
};

template <typename factor_t>
class huff_factor_iterator {
    bit_reader* r = nullptr;
    huffman* len_huff = nullptr;
    huffman* dist_huff = nullptr;
    uint64_t pos = 0;
    uint64_t n = 0;
    uint64_t remaining = 0;
    factor_t cur {};

    void advance()
    {
        if (pos >= n) return;

        if (remaining == 0) {
            remaining = get_elias_delta(*r);
            len_huff->read_table(*r, huff_sigma);
            dist_huff->read_table(*r, huff_sigma);
        }

        uint8_t lb = (uint8_t) len_huff->decode(*r);
        if (lb == 0) {
            cur.src = (decltype(cur.src)) r->get(8);
            cur.len = 0;
            pos++;
        } else {
            uint64_t len = ((uint64_t) 1 << (lb - 1)) | r->get(lb - 1);
            uint8_t db = (uint8_t) dist_huff->decode(*r);
            uint64_t dist = ((uint64_t) 1 << (db - 1)) | r->get(db - 1);
            cur.len = (decltype(cur.len)) len;
            cur.src = (decltype(cur.src)) (pos - dist);
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
        : r(&reader), len_huff(&lh), dist_huff(&dh), n(n) { advance(); }

    reference operator*() const { return cur; }
    huff_factor_iterator& operator++() { advance(); return *this; }
    huff_factor_iterator operator++(int) { auto tmp = *this; advance(); return tmp; }
};
