#pragma once

#include <cstdint>
#include <vector>

#include <tsl/sparse_set.h>
#include <lz77_sss/data_structures/rk61_substring.hpp>
#include <lce/lce_naive_wordwise_xor.hpp>
#include <lz77_sss/algorithms/lce_l.hpp>
#include <lz77_sss/misc/utils.hpp>

template <
    typename pos_t = uint32_t,
    typename sidx_t = uint32_t,
    typename lce_r_t = alx::lce::lce_naive_wordwise_xor<char>
> class sample_index {
    public:

    struct sxa_interval_t {
        sidx_t b;
        sidx_t e;
    };

    struct query_context {
        sidx_t b;
        sidx_t e;
        pos_t lce_b;
        pos_t lce_e;

        inline pos_t match_length() const {
            return std::min<pos_t>(lce_b, lce_e);
        }

        inline sxa_interval_t interval() const {
            return {b, e};
        }
    };

    protected:

    static sxa_interval_t pos_to_interval(pos_t pos) {
        sxa_interval_t iv;
        *reinterpret_cast<uint64_t*>(&iv) = uint64_t{pos} | (uint64_t{1} << 63);
        return iv;
    }

    template <direction dir>
    pos_t pos_from_interval(const sxa_interval_t& iv) const {
        const uint64_t& iv_u64 = *reinterpret_cast<const uint64_t*>(&iv);
        pos_t pos;

        if (iv_u64 & (uint64_t{1} << 63)) {
            pos = iv_u64 & (std::numeric_limits<uint64_t>::max() >> 1);
        } else {
            pos = sample(sxa<dir>(iv.b));
        }

        return pos;
    }

    static constexpr sidx_t no_occ = std::numeric_limits<sidx_t>::max();

    template <direction dir>
    struct sxivx_hash {
        sample_index* idx = nullptr;
        pos_t len, s;
        sxivx_hash() = default;
        sxivx_hash(sample_index* idx, pos_t len)
          : idx(idx), len(len), s(idx->num_samples()) {}

        inline std::size_t operator()(const sxa_interval_t& iv) const {
            return idx->rks.substring<dir>(idx->pos_from_interval<dir>(iv), len);
        }
    };

    template <direction dir>
    struct sxivx_eq {
        sample_index* idx = nullptr;
        pos_t len, s;
        sxivx_eq() = default;
        sxivx_eq(sample_index* idx, const pos_t len)
          : idx(idx), len(len), s(idx->num_samples()) {}

        inline bool operator()(const sxa_interval_t& iv_1, const sxa_interval_t& iv_2) const {
            pos_t pos_1 = idx->pos_from_interval<dir>(iv_1);
            pos_t pos_2 = idx->pos_from_interval<dir>(iv_2);
            return pos_1 == pos_2 || idx->lce<dir>(pos_1, pos_2, len) >= len;
        }
    };

    template <direction dir>
    using sxivx_map_t = tsl::sparse_set<
        sxa_interval_t,
        sxivx_hash<dir>,
        sxivx_eq<dir>
    >;

    pos_t n = 0;
    sidx_t s = 0;
    pos_t max_lce_l = 0;
    uint64_t byte_size = 0;

    const char* T = nullptr;
    const pos_t* S = nullptr;
    const lce_r_t* LCE_R = nullptr;
    rk61_substring rks;

    std::vector<sidx_t> SPA;
    std::vector<sidx_t> SSA;

    sxa_interval_t SCIV[1 << 8];
    sxa_interval_t SXIV2[2][1 << 16];
    std::vector<sxivx_map_t<LEFT>> SPIVX;
    std::vector<sxivx_map_t<RIGHT>> SSIVX;

    std::vector<pos_t> smpl_pat_lens[2];

    template <direction dir>
    inline std::vector<sxivx_map_t<dir>>& SXIVX() {
        if constexpr (dir == LEFT) {
            return SPIVX;
        } else {
            return SSIVX;
        }
    }

    template <direction dir>
    inline const std::vector<sxivx_map_t<dir>>& SXIVX() const {
        if constexpr (dir == LEFT) {
            return SPIVX;
        } else {
            return SSIVX;
        }
    }

    template <direction dir>
    inline sidx_t SXA(sidx_t i) const {
        if constexpr (dir == LEFT) {
            return SPA[i];
        } else {
            return SSA[i];
        }
    }

    template<direction dir>
    bool is_pos_in_T(pos_t p, pos_t offs) const {
        if constexpr (dir == LEFT) {
            return p >= offs;
        } else {
            return p + offs < n;
        }
    }

    template<typename T1, direction dir, pos_t offs, typename T2>
    static T1 val_offs(const T2* ptr, pos_t p) {
        if constexpr (dir == LEFT) {
            static constexpr pos_t offs_l = (sizeof(T1) - 1) + offs;
            return *reinterpret_cast<const T1*>(
                &ptr[p - offs_l]);
        } else {
            return *reinterpret_cast<const T1*>(
                &ptr[p + offs]);
        }
    }

    inline const sxa_interval_t& sciv(pos_t pos_char) const {
        return SCIV[val_offs<uint8_t, RIGHT, 0>(T, pos_char)];
    }

    template <direction dir>
    inline const sxa_interval_t& sxiv2(pos_t pos_pat) const {
        return SXIV2[dir][val_offs<uint16_t, dir, 0>(T, pos_pat)];
    }

    inline bool occurs1(pos_t pos_char) const {
        return sciv(pos_char).b != no_occ;
    }

    template <direction dir>
    inline bool occurs2(pos_t pos_pat) const {
        if constexpr (dir == LEFT) {
            if (pos_pat < 1) [[unlikely]] return false;
        } else {
            if (pos_pat >= n - 1) [[unlikely]] return false;
        }

        return sxiv2<dir>(pos_pat).b != no_occ;
    }

    template <direction dir>
    inline pos_t lce(
        pos_t i, pos_t j,
        pos_t max_lce_l = std::numeric_limits<pos_t>::max()
    ) const {
        if constexpr (dir == LEFT) {
            return lce_l_128<pos_t>(T, i, j, max_lce_l);
        } else {
            return LCE_R->lce(i, j);
        }
    }

    template <direction dir>
    inline pos_t lce_offs(
        pos_t i, pos_t j, pos_t offs,
        pos_t max_lce_l = std::numeric_limits<pos_t>::max()
    ) const {
        if constexpr (dir == LEFT) {
            if (!is_pos_in_T<dir>(std::min<pos_t>(i, j), offs)) [[unlikely]] {
                return offs;
            }

            return offs + lce_l_128<pos_t>(T, i - offs, j - offs, max_lce_l - offs);
        } else {
            if (!is_pos_in_T<dir>(std::max<pos_t>(i, j), offs)) [[unlikely]] {
                return offs;
            }

            return offs + LCE_R->lce(offs + i, offs + j);
        }
    }

    template <direction dir>
    inline bool cmp_lex(pos_t i, pos_t j, pos_t lce) const {
        if (i == j) [[unlikely]] {
            return false;
        }

        if constexpr (dir == LEFT) {
            if (lce > std::min<pos_t>(i, j)) [[unlikely]] {
                return i < j;
            }

            return char_to_uchar(T[i - lce]) < char_to_uchar(T[j - lce]);
        } else {
            if (std::max<pos_t>(i, j) + lce == n) [[unlikely]] {
                return i > j;
            }

            return char_to_uchar(T[i + lce]) < char_to_uchar(T[j + lce]);
        }
    }

    template <direction dir>
    inline bool cmp_sample_lex(sidx_t i, sidx_t j) const {
        if (i == j) [[unlikely]] {
            return false;
        }

        return cmp_lex<dir>(S[i], S[j], lce<dir>(S[i], S[j], max_lce_l));
    }

    void build_sxa12_intervals(uint16_t p, bool log);

    template <direction dir>
    void build_samples(pos_t typ_lce_r, uint16_t p, bool log);

    public:

    sample_index() = default;

    void build(
        const std::string& T,
        const std::vector<pos_t>& S,
        const lce_r_t& LCE_R,
        bool use_samples = true,
        uint16_t p = 1,
        bool log = false,
        pos_t max_lce_l = std::numeric_limits<pos_t>::max(),
        pos_t typ_lce_r = std::numeric_limits<pos_t>::max()
    ) {
        uint64_t baseline_memory_alloc = malloc_count_current();
        auto time = now();

        this->T = T.data();
        this->S = S.data();
        this->LCE_R = &LCE_R;
        this->max_lce_l = max_lce_l;
        n = T.size();
        s = S.size();

        #ifndef NDEBUG
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 1; i < s; i++) {
            assert(S[i] > S[i - 1]);
        }
        #endif

        if (log) {
            std::cout << "building SPA" << std::flush;
        }

        no_init_resize(SPA, s);

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < s; i++) {
            SPA[i] = i;
        }

        ips4o::parallel::sort(SPA.begin(), SPA.end(),
            [&](sidx_t i, sidx_t j){
                return cmp_sample_lex<LEFT>(i, j);
        });

        #ifndef NDEBUG
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 1; i < s; i++) {
            assert(!cmp_sample_lex<LEFT>(SPA[i], SPA[i - 1]));
        }
        #endif

        if (log) {
            std::cout << " (" << format_size(s * sizeof(sidx_t)) << ")" << std::flush;
            time = log_runtime(time);
            std::cout << "building SSA" << std::flush;
        }

        no_init_resize(SSA, s);
    
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < s; i++) {
            SSA[i] = i;
        }

        ips4o::parallel::sort(SSA.begin(), SSA.end(),
            [&](sidx_t i, sidx_t j){
                return cmp_sample_lex<RIGHT>(i, j);
        });

        #ifndef NDEBUG
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 1; i < s; i++) {
            assert(!cmp_sample_lex<RIGHT>(SSA[i], SSA[i - 1]));
        }
        #endif

        if (log) {
            std::cout << " (" << format_size(s * sizeof(sidx_t)) << ")" << std::flush;
            time = log_runtime(time);
        }

        if (use_samples) {
            if (log) {
                std::cout << "building rabin karp substring data structure" << std::flush;
            }

            uint64_t sampling_rate = max_lce_l < 256 ? n : max_lce_l;
            rks = rk61_substring(T, sampling_rate, 0, p);

            if (log) {
                std::cout << " (" << format_size(rks.size_in_bytes()) << ")" << std::flush;
                time = log_runtime(time);
            }

            build_sxa12_intervals(p, log);
            build_samples<LEFT>(typ_lce_r, p, log);
            build_samples<RIGHT>(typ_lce_r, p, log);
        }        

        byte_size = malloc_count_current() - baseline_memory_alloc;
    }

    const rk61_substring& rabin_karp_substring() const {
        return rks;
    }

    inline uint64_t size_in_bytes() const {
        return byte_size;
    }

    template <direction dir>
    inline sidx_t sxa(sidx_t i) const {
        if constexpr (dir == LEFT) {
            return SPA[i];
        } else {
            return SSA[i];
        }
    }

    inline sidx_t num_samples() const {
        return s;
    }

    inline pos_t sample(sidx_t i) const {
        return S[i];
    }

    inline sidx_t spa(sidx_t i) const {
        return SPA[i];
    }

    inline sidx_t ssa(sidx_t i) const {
        return SSA[i];
    }

    query_context query() const {
        return {
            .b = 0,
            .e = s - 1,
            .lce_b = 0,
            .lce_e = 0
        };
    }

    template <direction dir>
    inline query_context query(const sxa_interval_t& sxa_iv, pos_t pos_pat, pos_t len) const {
        return {
            .b = sxa_iv.b,
            .e = sxa_iv.e,
            .lce_b = lce_offs<dir>(pos_pat, S[SXA<dir>(sxa_iv.b)], len, max_lce_l),
            .lce_e = lce_offs<dir>(pos_pat, S[SXA<dir>(sxa_iv.e)], len, max_lce_l)
        };
    }

    inline query_context query_left(const sxa_interval_t& sxa_iv, pos_t pos_pat, pos_t len) const {
        return query<LEFT>(sxa_iv, pos_pat, len);
    }

    inline query_context query_right(const sxa_interval_t& sxa_iv, pos_t pos_pat, pos_t len) const {
        return query<RIGHT>(sxa_iv, pos_pat, len);
    }

    template <direction dir>
    inline pos_t num_sampled_pattern_lengths() const {
        return smpl_pat_lens[dir].size();
    }

    inline pos_t num_sampled_pattern_lengths_left() const {
        return smpl_pat_lens[LEFT].size();
    }

    inline pos_t num_sampled_pattern_lengths_right() const {
        return smpl_pat_lens[RIGHT].size();
    }

    template <direction dir>
    inline const std::vector<pos_t>& sampled_pattern_lengths() const {
        return smpl_pat_lens[dir];
    }

    inline const std::vector<pos_t>& sampled_pattern_lengths_left() const {
        return smpl_pat_lens[LEFT];
    }

    inline const std::vector<pos_t>& sampled_pattern_lengths_right() const {
        return smpl_pat_lens[RIGHT];
    }

    template <direction dir>
    inline std::pair<sxa_interval_t, bool> sxa_interval(
        pos_t pat_len_idx, pos_t pos_pat,
        std::size_t hash = std::numeric_limits<std::size_t>::max()
    ) const;

    inline std::pair<sxa_interval_t, bool> spa_interval(
        pos_t pat_len_idx, pos_t pos_pat,
        std::size_t hash = std::numeric_limits<std::size_t>::max()
    ) const {
        return sxa_interval<LEFT>(pat_len_idx, pos_pat, hash);
    }

    inline std::pair<sxa_interval_t, bool> ssa_interval(
        pos_t pat_len_idx, pos_t pos_pat,
        std::size_t hash = std::numeric_limits<std::size_t>::max()
    ) const {
        return sxa_interval<RIGHT>(pat_len_idx, pos_pat, hash);
    }

    template <direction dir>
    bool extend(const query_context& qc_old, query_context& qc_new, pos_t pos_pat, pos_t len, bool use_samples = true) const;
    
    bool extend_left(const query_context& qc_old, query_context& qc_new, pos_t pos_pat, pos_t len, bool use_samples = true) const {
        return extend<LEFT>(qc_old, qc_new, pos_pat, len, use_samples);
    }
    
    bool extend_right(const query_context& qc_old, query_context& qc_new, pos_t pos_pat, pos_t len, bool use_samples = true) const {
        return extend<RIGHT>(qc_old, qc_new, pos_pat, len, use_samples);
    }
    
    bool extend_left(query_context& qc, pos_t pos_pat, pos_t len, bool use_samples = true) const {
        return extend<LEFT>(qc, qc, pos_pat, len, use_samples);
    }
    
    bool extend_right(query_context& qc, pos_t pos_pat, pos_t len, bool use_samples = true) const {
        return extend<RIGHT>(qc, qc, pos_pat, len, use_samples);
    }

    template <direction dir>
    bool extend(query_context& qc, pos_t pos_pat, pos_t len, bool use_samples = true) const {
        return extend<dir>(qc, qc, pos_pat, len, use_samples);
    }

    template <direction dir>
    query_context interpolate(const query_context& qc_short, const query_context& qc_long, pos_t pos_pat, pos_t len) const;

    query_context interpolate_left(const query_context& qc_short, const query_context& qc_long, pos_t pos_pat, pos_t len) const {
        return interpolate<LEFT>(qc_short, qc_long, pos_pat, len);
    };

    query_context interpolate_right(const query_context& qc_short, const query_context& qc_long, pos_t pos_pat, pos_t len) const {
        return interpolate<RIGHT>(qc_short, qc_long, pos_pat, len);
    };

    template <direction dir>
    void locate(const query_context& qc, std::vector<pos_t>& Occ) const {
        Occ.reserve(Occ.size() + qc.e - qc.b + 1);

        for (sidx_t i = qc.b; i <= qc.e; i++) {
            Occ.emplace_back(S[SXA<dir>(i)]);
        }
    }

    template <direction dir>
    std::vector<pos_t> locate(const query_context& qc) const {
        std::vector<pos_t> Occ;
        locate(qc, Occ);
        return Occ;
    }
};

#include "construction.cpp"
#include "queries.cpp"