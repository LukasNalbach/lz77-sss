#pragma once

#include <cstdint>
#include <vector>

#include <ds/lce_naive_wordwise_xor.hpp>
#include <lz77_sss/algorithms/lce_l.hpp>
#include <lz77_sss/data_structures/rabin_karp_substring.hpp>
#include <lz77_sss/misc/utils.hpp>
#include <tsl/sparse_set.h>

template <
    typename pos_t = uint32_t,
    typename sidx_t = uint32_t,
    typename char_t = char,
    typename lce_r_t = lce::ds::lce_naive_wordwise_xor<char_t>>
class sample_index {
public:
    struct interval_t {
        sidx_t b;
        sidx_t e;

        bool empty() const
        {
            return b > e;
        }
    };

    struct query_ctx_t {
        sidx_t b;
        sidx_t e;
        pos_t lce_b;
        pos_t lce_e;

        inline pos_t match_length() const
        {
            return std::min<pos_t>(lce_b, lce_e);
        }

        inline interval_t interval() const
        {
            return { b, e };
        }
    };

protected:
    inline static interval_t pos_to_interval(pos_t pos)
    {
        interval_t iv;
        *reinterpret_cast<uint64_t*>(&iv) = uint64_t { pos } | (uint64_t { 1 } << 63);
        return iv;
    }

    template <direction dir>
    inline pos_t interval_to_pos(const interval_t& iv) const
    {
        const uint64_t& iv_u64 = *reinterpret_cast<const uint64_t*>(&iv);
        pos_t pos;

        if (iv_u64 & (uint64_t { 1 } << 63)) {
            pos = iv_u64 & (std::numeric_limits<uint64_t>::max() >> 1);
        } else {
            pos = sample(xa_s<dir>(iv.b));
        }

        return pos;
    }

    static constexpr sidx_t no_occ = std::numeric_limits<sidx_t>::max();

    template <direction dir>
    struct xiv_s_hash {
        sample_index* idx = nullptr;
        pos_t len;
        xiv_s_hash() = default;
        xiv_s_hash(sample_index* idx, pos_t len) : idx(idx), len(len) { }

        inline std::size_t operator()(const interval_t& iv) const
        {
            return idx->rks.template substring<dir>(idx->interval_to_pos<dir>(iv), len);
        }
    };

    template <direction dir>
    struct xiv_s_eq {
        sample_index* idx = nullptr;
        pos_t len;
        xiv_s_eq() = default;
        xiv_s_eq(sample_index* idx, const pos_t len) : idx(idx), len(len) { }

        inline bool operator()(const interval_t& iv_1, const interval_t& iv_2) const
        {
            pos_t pos_1 = idx->interval_to_pos<dir>(iv_1);
            pos_t pos_2 = idx->interval_to_pos<dir>(iv_2);
            return pos_1 == pos_2 || idx->lce<dir>(pos_1, pos_2, len) >= len;
        }
    };

    template <direction dir>
    using xiv_s_map_t = tsl::sparse_set<interval_t, xiv_s_hash<dir>, xiv_s_eq<dir>>;

    pos_t n = 0;
    sidx_t s = 0;
    pos_t max_patt_len_left = 0;
    uint64_t byte_size = 0;

    const char_t* T = nullptr;
    const pos_t* S = nullptr;
    const lce_r_t* LCE_R = nullptr;
    rabin_karp_substring<31, pos_t> rks;

    std::vector<sidx_t> PA_S;
    std::vector<sidx_t> SA_S;

    interval_t SIV_S_1[1 << 8];
    interval_t XIV_S_2[2][1 << 16];
    std::vector<xiv_s_map_t<LEFT>> PIV_S;
    std::vector<xiv_s_map_t<RIGHT>> SIV_S;

    std::vector<pos_t> smpl_pat_lens[2];

    template <direction dir>
    inline std::vector<xiv_s_map_t<dir>>& XIV_S()
    {
        if constexpr (dir == LEFT) {
            return PIV_S;
        } else {
            return SIV_S;
        }
    }

    template <direction dir>
    inline const std::vector<xiv_s_map_t<dir>>& XIV_S() const
    {
        if constexpr (dir == LEFT) {
            return PIV_S;
        } else {
            return SIV_S;
        }
    }

    template <direction dir>
    inline sidx_t XA_S(sidx_t i) const
    {
        if constexpr (dir == LEFT) {
            return PA_S[i];
        } else {
            return SA_S[i];
        }
    }

    template <direction dir>
    inline bool is_pos_in_T(pos_t p, pos_t offs) const
    {
        if constexpr (dir == LEFT) {
            return p >= offs;
        } else {
            return p + offs < n;
        }
    }

    template <typename T1, direction dir, typename T2>
    inline static T1 val_offs(const T2* ptr, pos_t p)
    {
        if constexpr (dir == LEFT) {
            return *reinterpret_cast<const T1*>(&ptr[p - (sizeof(T1) - 1)]);
        } else {
            return *reinterpret_cast<const T1*>(&ptr[p]);
        }
    }

    inline const interval_t& siv_s_1(pos_t pos_char) const
    {
        return SIV_S_1[val_offs<uint8_t, RIGHT>(T, pos_char)];
    }

    template <direction dir>
    inline const interval_t& xiv_s_2(pos_t pos_patt) const
    {
        return XIV_S_2[dir][val_offs<uint16_t, dir>(T, pos_patt)];
    }

    inline bool occurs1(pos_t pos_char) const
    {
        return siv_s_1(pos_char).b != no_occ;
    }

    template <direction dir>
    inline bool occurs2(pos_t pos_patt) const
    {
        if constexpr (dir == LEFT) {
            if (pos_patt < 1) [[unlikely]] return false;
        } else {
            if (pos_patt >= n - 1) [[unlikely]] return false;
        }

        return xiv_s_2<dir>(pos_patt).b != no_occ;
    }

    template <direction dir>
    inline pos_t lce(
        pos_t i, pos_t j,
        pos_t max_patt_len_left = std::numeric_limits<pos_t>::max()) const
    {
        if constexpr (dir == LEFT) {
            return lce_l_64<pos_t>(T, i, j, max_patt_len_left);
        } else {
            return LCE_R->lce(i, j);
        }
    }

    template <direction dir>
    inline pos_t lce_offs(
        pos_t i, pos_t j, pos_t offs,
        pos_t max_patt_len_left = std::numeric_limits<pos_t>::max()) const
    {
        if constexpr (dir == LEFT) {
            if (!is_pos_in_T<dir>(std::min<pos_t>(i, j), offs)) [[unlikely]] return offs;
            return offs + lce_l_64<pos_t>(T, i - offs, j - offs, max_patt_len_left - offs);
        } else {
            if (!is_pos_in_T<dir>(std::max<pos_t>(i, j), offs)) [[unlikely]] return offs;
            return offs + LCE_R->lce(offs + i, offs + j);
        }
    }

    template <direction dir>
    inline bool cmp_lex(pos_t i, pos_t j, pos_t lce) const
    {
        if (i == j) [[unlikely]] return false;

        if constexpr (dir == LEFT) {
            if (lce > std::min<pos_t>(i, j)) [[unlikely]] return i < j;
            return char_to_uchar(T[i - lce]) < char_to_uchar(T[j - lce]);
        } else {
            if (std::max<pos_t>(i, j) + lce == n) [[unlikely]] return i > j;
            return char_to_uchar(T[i + lce]) < char_to_uchar(T[j + lce]);
        }
    }

    template <direction dir>
    inline bool cmp_sample_lex(sidx_t i, sidx_t j) const
    {
        if (i == j) [[unlikely]] return false;
        return cmp_lex<dir>(S[i], S[j], lce<dir>(S[i], S[j], max_patt_len_left));
    }

    void build_xa_s_1_2_intervals(uint16_t p, bool log);

    template <direction dir>
    void build_samples(pos_t max_smpl_len, uint16_t p, bool log);

public:
    sample_index() = default;

    void build(
        const char_t* T,
        pos_t n,
        const std::vector<pos_t>& S,
        const lce_r_t& LCE_R,
        bool build_interval_samples = true,
        pos_t rks_sample_rate = 32,
        uint16_t p = 1,
        bool log = false,
        pos_t max_patt_len_left = std::numeric_limits<pos_t>::max(),
        pos_t max_smpl_len_right = std::numeric_limits<pos_t>::max())
    {
        uint64_t baseline_memory_alloc = malloc_count_current();
        auto time = now();

        this->T = T;
        this->S = S.data();
        this->LCE_R = &LCE_R;
        this->max_patt_len_left = max_patt_len_left;
        this->n = n;
        s = S.size();

        #ifndef NDEBUG
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 1; i < s; i++) {
            assert(S[i] > S[i - 1]);
        }
        #endif

        if (log) {
            std::cout << "building PA_C" << std::flush;
        }

        no_init_resize(PA_S, s);

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < s; i++) {
            PA_S[i] = i;
        }

        ips4o::parallel::sort(PA_S.begin(), PA_S.end(),
            [&](sidx_t i, sidx_t j) {
                return cmp_sample_lex<LEFT>(i, j);
            });

        #ifndef NDEBUG
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 1; i < s; i++) {
            assert(!cmp_sample_lex<LEFT>(PA_S[i], PA_S[i - 1]));
        }
        #endif

        if (log) {
            log_phase("pa_c", time_diff_ns(time, now()));
            std::cout << " (" << format_size(s * sizeof(sidx_t)) << ")" << std::flush;
            time = log_runtime(time);
            std::cout << "building SA_C" << std::flush;
        }

        no_init_resize(SA_S, s);

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < s; i++) {
            SA_S[i] = i;
        }

        ips4o::parallel::sort(SA_S.begin(), SA_S.end(),
            [&](sidx_t i, sidx_t j) {
                return cmp_sample_lex<RIGHT>(i, j);
            });

        #ifndef NDEBUG
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 1; i < s; i++) {
            assert(!cmp_sample_lex<RIGHT>(SA_S[i], SA_S[i - 1]));
        }
        #endif

        if (log) {
            log_phase("sa_c", time_diff_ns(time, now()));
            std::cout << " (" << format_size(s * sizeof(sidx_t)) << ")" << std::flush;
            time = log_runtime(time);
        }

        if (build_interval_samples) {
            if (log) {
                std::cout << "building rabin karp substring data structure" << std::flush;
            }

            rks = rabin_karp_substring<31, pos_t>(T, n, rks_sample_rate, 0, p);

            if (log) {
                log_phase("rks", time_diff_ns(time, now()));
                std::cout << " (" << format_size(rks.size_in_bytes()) << ")" << std::flush;
                time = log_runtime(time);
            }

            build_xa_s_1_2_intervals(p, log);
            build_samples<LEFT>(max_patt_len_left, p, log);
            build_samples<RIGHT>(max_smpl_len_right, p, log);
        }

        byte_size = malloc_count_current() - baseline_memory_alloc;
    }

    inline const rabin_karp_substring<31, pos_t>& rab_karp_substr() const
    {
        return rks;
    }

    inline uint64_t size_in_bytes() const
    {
        return byte_size;
    }

    template <direction dir>
    inline sidx_t xa_s(sidx_t i) const
    {
        if constexpr (dir == LEFT) {
            return PA_S[i];
        } else {
            return SA_S[i];
        }
    }

    inline sidx_t num_samples() const
    {
        return s;
    }

    inline pos_t sample(sidx_t i) const
    {
        return S[i];
    }

    inline sidx_t pa_s(sidx_t i) const
    {
        return PA_S[i];
    }

    inline sidx_t sa_s(sidx_t i) const
    {
        return SA_S[i];
    }

    inline query_ctx_t query() const
    {
        return {
            .b = 0,
            .e = s - 1,
            .lce_b = 0,
            .lce_e = 0
        };
    }

    template <direction dir>
    inline query_ctx_t query(const interval_t& sxa_iv, pos_t pos_patt, pos_t len) const
    {
        return {
            .b = sxa_iv.b,
            .e = sxa_iv.e,
            .lce_b = lce_offs<dir>(pos_patt, S[XA_S<dir>(sxa_iv.b)], len, max_patt_len_left),
            .lce_e = lce_offs<dir>(pos_patt, S[XA_S<dir>(sxa_iv.e)], len, max_patt_len_left)
        };
    }

    inline query_ctx_t query_left(const interval_t& sxa_iv, pos_t pos_patt, pos_t len) const
    {
        return query<LEFT>(sxa_iv, pos_patt, len);
    }

    inline query_ctx_t query_right(const interval_t& sxa_iv, pos_t pos_patt, pos_t len) const
    {
        return query<RIGHT>(sxa_iv, pos_patt, len);
    }

    template <direction dir>
    inline pos_t num_sampled_pattern_lengths() const
    {
        return smpl_pat_lens[dir].size();
    }

    inline pos_t num_sampled_pattern_lengths_left() const
    {
        return smpl_pat_lens[LEFT].size();
    }

    inline pos_t num_sampled_pattern_lengths_right() const
    {
        return smpl_pat_lens[RIGHT].size();
    }

    template <direction dir>
    inline const std::vector<pos_t>& sampled_pattern_lengths() const
    {
        return smpl_pat_lens[dir];
    }

    inline const std::vector<pos_t>& sampled_pattern_lengths_left() const
    {
        return smpl_pat_lens[LEFT];
    }

    inline const std::vector<pos_t>& sampled_pattern_lengths_right() const
    {
        return smpl_pat_lens[RIGHT];
    }

    template <direction dir>
    inline std::pair<interval_t, bool> sxa_interval(
        pos_t pat_len_idx, pos_t pos_patt,
        std::size_t hash = std::numeric_limits<std::size_t>::max()) const;

    inline std::pair<interval_t, bool> pa_s_interval(
        pos_t pat_len_idx, pos_t pos_patt,
        std::size_t hash = std::numeric_limits<std::size_t>::max()) const
    {
        return sxa_interval<LEFT>(pat_len_idx, pos_patt, hash);
    }

    inline std::pair<interval_t, bool> sa_s_interval(
        pos_t pat_len_idx, pos_t pos_patt,
        std::size_t hash = std::numeric_limits<std::size_t>::max()) const
    {
        return sxa_interval<RIGHT>(pat_len_idx, pos_patt, hash);
    }

    template <direction dir>
    inline bool extend(const query_ctx_t& qc_old, query_ctx_t& qc_new, pos_t pos_patt, pos_t len, bool use_interval_samples = true) const;

    bool extend_left(const query_ctx_t& qc_old, query_ctx_t& qc_new, pos_t pos_patt, pos_t len, bool use_interval_samples = true) const
    {
        return extend<LEFT>(qc_old, qc_new, pos_patt, len, use_interval_samples);
    }

    inline bool extend_right(const query_ctx_t& qc_old, query_ctx_t& qc_new, pos_t pos_patt, pos_t len, bool use_interval_samples = true) const
    {
        return extend<RIGHT>(qc_old, qc_new, pos_patt, len, use_interval_samples);
    }

    inline bool extend_left(query_ctx_t& qc, pos_t pos_patt, pos_t len, bool use_interval_samples = true) const
    {
        return extend<LEFT>(qc, qc, pos_patt, len, use_interval_samples);
    }

    inline bool extend_right(query_ctx_t& qc, pos_t pos_patt, pos_t len, bool use_interval_samples = true) const
    {
        return extend<RIGHT>(qc, qc, pos_patt, len, use_interval_samples);
    }

    template <direction dir>
    inline bool extend(query_ctx_t& qc, pos_t pos_patt, pos_t len, bool use_interval_samples = true) const
    {
        return extend<dir>(qc, qc, pos_patt, len, use_interval_samples);
    }

    template <direction dir>
    query_ctx_t interpolate(const query_ctx_t& qc_short, const query_ctx_t& qc_long, pos_t pos_patt, pos_t len) const;

    inline query_ctx_t interpolate_left(const query_ctx_t& qc_short, const query_ctx_t& qc_long, pos_t pos_patt, pos_t len) const
    {
        return interpolate<LEFT>(qc_short, qc_long, pos_patt, len);
    };

    inline query_ctx_t interpolate_right(const query_ctx_t& qc_short, const query_ctx_t& qc_long, pos_t pos_patt, pos_t len) const
    {
        return interpolate<RIGHT>(qc_short, qc_long, pos_patt, len);
    };

    template <direction dir>
    inline void locate(const query_ctx_t& qc, std::vector<pos_t>& Occ) const
    {
        Occ.reserve(Occ.size() + qc.e - qc.b + 1);

        for (sidx_t i = qc.b; i <= qc.e; i++) {
            Occ.emplace_back(S[XA_S<dir>(i)]);
        }
    }

    template <direction dir>
    inline std::vector<pos_t> locate(const query_ctx_t& qc) const
    {
        std::vector<pos_t> Occ;
        locate(qc, Occ);
        return Occ;
    }
};

#include "construction.cpp"
#include "queries.cpp"