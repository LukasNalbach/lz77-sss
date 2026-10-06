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

#include <filesystem>
#include <functional>
#include <vector>

#include <ds/lce_sss.hpp>
#include <text/direct_text.hpp>
#include <text/packed_text.hpp>
#include <text/split_text.hpp>

#include <lz77_sss/data_structures/range/range.hpp>
#include <lz77_sss/data_structures/min_tree.hpp>
#include <lz77_sss/data_structures/exact_gap_index.hpp>
#include <lz77_sss/data_structures/rolling_hash_index.hpp>
#include <lz77_sss/data_structures/sample_index/sample_index.hpp>
#include <lz77_sss/misc/direct_io.hpp>
#include <lz77_sss/data_structures/bit_aligned_interleaved_vectors.hpp>
#include <lz77_sss/misc/progress.hpp>
#include <lz77_sss/misc/log.hpp>
#include <lz77_sss/misc/search.hpp>
#include <lz77_sss/misc/utils.hpp>

class lz77_sss {
public:
    enum factorize_mode {
        auto_gaps = 1,
        skip_gaps = 2,
        exact_gaps = 3
    };

    enum transform_mode {
        with_interval_samples = 1,
        without_interval_samples = 2
    };

    enum exact_algorithm {
        auto_select = 1,
        sss_based = 2,
        sa_based = 3
    };

    struct parameters {
        uint16_t num_threads = 0;
        bool log = false;
        std::chrono::steady_clock::time_point time_start { };
        uint64_t tau = 512;
        factorize_mode fact_mode = auto_gaps;
        transform_mode transf_mode = with_interval_samples;
        range_ds_kind range_ds { .type = range_ds_type::swsg, .decomposed = true };
        exact_algorithm exact_alg = auto_select;
        bool log_factor_count = true;
    };

    static constexpr uint64_t       max_input_size                      = uint40_t::max;

    static constexpr factorize_mode default_fact_mode                   = auto_gaps;
    static constexpr transform_mode default_transf_mode                 = with_interval_samples;

    static constexpr uint64_t       default_tau                         = 512;
    static constexpr uint64_t       max_delta                           = 256;
    static constexpr uint64_t       rks_sample_rate                     = 32;
    static constexpr uint64_t       range_scan_threshold                = 4096;
    static constexpr uint64_t       par_gap_blk_len                     = 4096;
    static constexpr uint64_t       par_blks_per_thr                    = 4;
    static constexpr uint64_t       par_owners_per_thr                  = 4;
    static constexpr uint64_t       par_fact_bufs                       = 4;
    static constexpr uint64_t       sa_sects_ahead_per_thr              = 4;
    static constexpr uint64_t       num_patt_lens                       = 5;
    static constexpr uint64_t       par_sects_per_thr                   = 16;
    static constexpr uint64_t       lpf_chunks_per_thr                  = 64;
    static constexpr uint64_t       max_sa_input_size                   = 1ULL << 38;
    static constexpr uint64_t       min_sa_sect_len                     = 64;
    static constexpr uint64_t       max_sa_sect_len                     = 1'048'576;
    static constexpr uint64_t       max_gap_ctx_per_tau                 = 2;
    static constexpr uint64_t       min_hash_gap_ctx                    = 16;
    static constexpr double         hash_gap_ctx_factor                 = 24;
    static constexpr double         max_sa_peak_ratio                   = 1.25;
    static constexpr double         sa_extra_bytes_per_char             = 0.5;
    static constexpr double         without_iv_smpl_peak_bytes_per_char = 1.7;
    static constexpr double         without_iv_smpl_peak_bytes_per_fact = 44.6;
    static constexpr double         with_iv_smpl_peak_bytes_per_char    = 1.34;
    static constexpr double         with_iv_smpl_peak_bytes_per_fact    = 115.7;

    using patt_len_entry_t = std::pair<double, std::array<uint64_t, num_patt_lens>>;
    static constexpr double         infty                               = std::numeric_limits<double>::max();

    static constexpr std::array<patt_len_entry_t, 10> patt_len_table {
        patt_len_entry_t { 6, { 2, 3, 4, 5, 6 } },
        patt_len_entry_t { 8, { 2, 3, 4, 6, 8 } },
        patt_len_entry_t { 12, { 2, 3, 4, 8, 12 } },
        patt_len_entry_t { 16, { 2, 4, 6, 9, 16 } },
        patt_len_entry_t { 32, { 2, 4, 6, 10, 20 } },
        patt_len_entry_t { 64, { 2, 4, 7, 12, 28 } },
        patt_len_entry_t { 128, { 2, 4, 8, 16, 36 } },
        patt_len_entry_t { 256, { 2, 5, 10, 20, 42 } },
        patt_len_entry_t { 1024, { 2, 6, 12, 24, 48 } },
        patt_len_entry_t { infty, { 2, 8, 16, 32, 64 } }
    };
    
    static double get_patt_len_guess(double avg_gap_len, double avg_lpf_phr_len, double rel_len_gaps)
    {
        return std::min<double>({avg_gap_len, avg_lpf_phr_len, 8.0 * std::pow(128, 1.0 - rel_len_gaps)});
    }

    static uint64_t exact_gaps_bytes(uint64_t len_G, uint64_t num_gaps, uint64_t symbol_bytes = 1)
    {
        const uint64_t sa_bytes = len_G <= INT32_MAX ? sizeof(int32_t) : sizeof(sa_int40_t);
        const uint64_t gap_bytes = symbol_bytes == 1 ? 1 : sa_bytes;
        return len_G * (gap_bytes + 2 * sa_bytes) + (len_G / 8) * (1 + sizeof(factor)) +
            num_gaps * 3 * sizeof(uint64_t);
    }

    static double estimated_sa_peak(uint64_t n, uint64_t sa_bytes, uint64_t other_bytes)
    {
        return other_bytes + n * (2.0 * sa_bytes + sa_extra_bytes_per_char);
    }

    static double estimated_sss_peak(uint64_t n, uint64_t num_aprx_fact, transform_mode transf_mode)
    {
        return transf_mode == with_interval_samples
            ? n * with_iv_smpl_peak_bytes_per_char + num_aprx_fact * with_iv_smpl_peak_bytes_per_fact
            : n * without_iv_smpl_peak_bytes_per_char + num_aprx_fact * without_iv_smpl_peak_bytes_per_fact;
    }

    static double estimated_sss_peak(uint64_t n, uint64_t text_bytes, uint64_t num_aprx_fact, transform_mode transf_mode)
    {
        return estimated_sss_peak(n, num_aprx_fact, transf_mode) + double(text_bytes) - double(n);
    }

    static uint64_t get_max_smpl_len_right(double aprx_comp_ratio)
    {
        return std::round(aprx_comp_ratio * (1.0 + 0.5 * std::exp(-aprx_comp_ratio / 1000.0)));
    }

    struct factor {
        uint64_t src;
        uint64_t len;

        friend class lz77_sss;

        static factor literal(uint64_t chr) { return factor { .src = chr, .len = 0 }; }

        static factor gap(uint64_t len) { return factor { .src = len, .len = 0 }; }

        bool is_literal() const { return len == 0; }

        bool is_gap() const { return len == 0; }

        uint64_t text_len() const { return std::max<uint64_t>(1, len); }

        static uint8_t bytes_for(uint64_t max_value) { return std::max<uint8_t>(1, (std::bit_width(max_value) + 7) / 8); }

        std::ostream& write(std::ostream& out, uint8_t src_bytes, uint8_t len_bytes) const
        {
            out.write((const char*) &src, src_bytes);
            return out.write((const char*) &len, len_bytes);
        }

        std::istream& read(std::istream& in, uint8_t src_bytes, uint8_t len_bytes)
        {
            src = 0;
            len = 0;
            in.read((char*) &src, src_bytes);
            return in.read((char*) &len, len_bytes);
        }
    };

    class factor_sink {
        static constexpr uint64_t capacity = 8192;

        void* ctx = nullptr;
        void (*emit)(void*, const factor*, uint64_t) = nullptr;
        std::vector<factor> buf;

    public:
        factor_sink() = default;

        template <typename fnc_t>
        explicit factor_sink(fnc_t& fnc)
            : ctx(&fnc)
            , emit([](void* c, const factor* f, uint64_t num) {
                  fnc_t& g = *static_cast<fnc_t*>(c);
                  for (uint64_t i = 0; i < num; i++) g(f[i]);
              })
        {
            buf.reserve(capacity);
        }

        inline void operator()(factor f)
        {
            buf.emplace_back(f);
            if (buf.size() == capacity) [[unlikely]] flush();
        }

        void flush()
        {
            if (buf.empty()) return;
            emit(ctx, buf.data(), buf.size());
            buf.clear();
        }
    };

    using direct_text = lce::text::direct_text<char>;
    using packed_text = lce::text::packed_text<>;
    using split_text = lce::text::split_text<>;
    using int16_direct_text = lce::text::direct_text<uint16_t>;
    using int_direct_text = lce::text::direct_text<uint32_t>;
    using int40_direct_text = lce::text::direct_text<uint40_t>;
    using int64_direct_text = lce::text::direct_text<uint64_t>;
    using int_packed_text = lce::text::packed_text<uint32_t>;
    using int_split_text = lce::text::split_text<uint32_t>;

    template <typename text_t, typename output_fnc_t>
    static void factorize_approximate(text_t text, output_fnc_t output, parameters params = { })
    {
        factor_sink sink(output);
        factorizer<text_t>(text, params).factorize(approximate, sink);
        sink.flush();
    }

    template <typename text_t, typename output_fnc_t>
    static void factorize_exact(text_t text, output_fnc_t output, parameters params = { })
    {
        factor_sink sink(output);
        factorizer<text_t>(text, params).factorize(exact, sink);
        sink.flush();
    }

    template <typename text_t, typename gapped_fnc_t>
    static void factorize_gapped(text_t text, gapped_fnc_t fnc, parameters params = { })
    {
        params.fact_mode = skip_gaps;
        factorizer<text_t> impl(text, params);
        impl.gapped_fnc = [&fnc](const gapped_factorization& gapped) { fnc(gapped); };
        factor_sink sink;
        impl.factorize(approximate, sink);
    }

    template <typename output_fnc_t>
    static void factorize_approximate(char* input, uint64_t input_size, output_fnc_t output, parameters params = { })
    {
        factorize_approximate(direct_text(input, input_size), output, params);
    }

    template <typename output_fnc_t>
    static void factorize_exact(char* input, uint64_t input_size, output_fnc_t output, parameters params = { })
    {
        factorize_exact(direct_text(input, input_size), output, params);
    }

    template <typename int_t>
    static constexpr bool is_int_input_v = std::is_same_v<int_t, uint40_t> || (std::is_integral_v<int_t> &&
        !std::is_const_v<int_t> && !std::is_same_v<int_t, bool> &&
        (sizeof(int_t) == 1 || sizeof(int_t) == 2 || sizeof(int_t) == 4 || sizeof(int_t) == 8));

    template <typename int_t>
    using int_symbol_t = typename std::conditional_t<std::is_same_v<int_t, uint40_t>, std::type_identity<uint40_t>,
        std::conditional_t<sizeof(int_t) == 1, std::type_identity<char>, std::make_unsigned<int_t>>>::type;

    template <typename int_t>
    using int_text_t = lce::text::direct_text<int_symbol_t<int_t>>;

    template <typename int_t, typename output_fnc_t>
    static void factorize_approximate(int_t* input, uint64_t input_size, output_fnc_t output, parameters params = { })
        requires(is_int_input_v<int_t>);

    template <typename int_t, typename output_fnc_t>
    static void factorize_exact(int_t* input, uint64_t input_size, output_fnc_t output, parameters params = { })
        requires(is_int_input_v<int_t>);

    template <typename int_t>
    static int_text_t<int_t> int_text(int_t* input, uint64_t input_size, uint16_t p = 0)
        requires(is_int_input_v<int_t>);
    
    template <std::input_iterator fact_it_t, typename out_it_t>
    static void decode(fact_it_t fact_it, out_it_t out_it, uint64_t output_size);

protected:
    struct __attribute__((packed)) lpf_phrase {
        uint40_t beg;
        uint40_t end;
        uint40_t src;
    };

    class lpf_array {
    public:
        void init(uint64_t max_value) { data.reset(0, { max_value, max_value, max_value }); }

        inline void emplace_back(lpf_phrase phrase)
        {
            data.push_back({ uint64_t(phrase.beg), uint64_t(phrase.end), uint64_t(phrase.src) });
        }

        inline lpf_phrase operator[](uint64_t i) const
        {
            return lpf_phrase { .beg = uint40_t(data.get<0>(i)), .end = uint40_t(data.get<1>(i)),
                .src = uint40_t(data.get<2>(i)) };
        }

        inline void set(uint64_t i, lpf_phrase phrase)
        {
            data.set<0>(i, uint64_t(phrase.beg));
            data.set<1>(i, uint64_t(phrase.end));
            data.set<2>(i, uint64_t(phrase.src));
        }

        inline lpf_phrase back() const { return (*this)[data.size() - 1]; }
        inline uint64_t size() const { return data.size(); }
        inline bool empty() const { return data.size() == 0; }
        uint64_t size_in_bytes() const { return data.size_in_bytes(); }

        std::vector<lpf_phrase> to_vector() const
        {
            std::vector<lpf_phrase> out;
            no_init_resize(out, data.size());
            for (uint64_t i = 0; i < data.size(); i++) out[i] = (*this)[i];
            return out;
        }

        void from_vector(std::vector<lpf_phrase>& in, uint64_t max_value)
        {
            data.reset(in.size(), { max_value, max_value, max_value });
            for (const lpf_phrase& phrase : in) emplace_back(phrase);
            in.clear();
            in.shrink_to_fit();
        }

    private:
        bit_aligned_interleaved_vectors<3> data;
    };

public:
    class gapped_factorization {
    public:
        gapped_factorization(std::vector<lpf_array>& phrases, uint64_t n);

        uint64_t num_sections() const { return sections.size(); }

        uint64_t section_begin(uint64_t s) const { return s == 0 ? 0 : sections[s].beg; }

        void emit_section(uint64_t s, factor_sink& out) const;

    private:
        struct sect_t {
            uint64_t chunk;
            uint64_t first;
            uint64_t beg;
            uint64_t nxt;
        };

        const std::vector<lpf_array>* lpf_chunks = nullptr;
        uint64_t n = 0;
        std::vector<sect_t> sections;
    };

protected:
    enum quality_mode {
        approximate,
        exact
    };

    lz77_sss() = delete;
    lz77_sss(lz77_sss&& other) = delete;
    lz77_sss(const lz77_sss& other) = delete;
    lz77_sss& operator=(lz77_sss&& other) = delete;
    lz77_sss& operator=(const lz77_sss& other) = delete;

    template <typename text_t>
    class factorizer {
    public:
        using lce_t = lce::ds::lce_sss<text_t, uint40_t>;
        using rh_idx_t = rolling_hash_index<num_patt_lens, text_t>;
        using symbol_t = typename text_t::symbol_type;

        static constexpr uint64_t symbol_bytes = text_t::is_byte_text ? 1 : sizeof(uint32_t);

        static uint8_t src_bytes_of(const text_t& T)
        {
            return factor::bytes_for(std::max<uint64_t>(T.size(), text_t::is_byte_text ? 256 : T.sigma()) - 1);
        }

        static uint8_t len_bytes_of(const text_t& T) { return factor::bytes_for(T.size()); }

        static uint64_t target_rh_idx_bytes_for(uint64_t n, double rel_len_gaps)
        {
            return std::min<uint64_t>(rh_idx_t::max_bytes, std::max<uint64_t>({rh_idx_t::min_bytes,
                malloc_count_peak() - malloc_count_current(),
                (uint64_t)((n / 3.0) * rel_len_gaps)}));
        }

        std::chrono::steady_clock::time_point time_start, time;
        uint64_t baseline_bytes = 0;
        uint64_t target_rh_idx_bytes = 0;
        bool log = false;
        uint16_t p = 0;
        uint64_t tau = default_tau;
        factorize_mode fact_mode = default_fact_mode;
        transform_mode transf_mode = default_transf_mode;
        range_ds_kind kind { .type = range_ds_type::swsg, .decomposed = true };
        exact_algorithm exact_alg = auto_select;
        bool log_factor_count = true;
        std::function<void(const gapped_factorization&)> gapped_fnc;

        text_t T;
        uint64_t n = 0;
        uint64_t size_sss = 0;
        uint64_t num_lpf_phr = 0;
        uint64_t len_lpf_phr = 0;
        uint64_t num_fact = 0;
        uint64_t len_gaps = 0;
        uint64_t num_gaps = 0;
        uint64_t gap_ctx = 0;
        bool factorize_gaps_exact = false;

        lce_t LCE;
        std::vector<lpf_array> LPF;
        std::array<uint64_t, num_patt_lens> patt_lens;
        rh_idx_t rh_idx;


        struct lpf_cursor_t {
            uint64_t chunk;
            uint64_t i;
        };
        
        struct par_blk_t {
            uint64_t beg;
            lpf_cursor_t lpf_cursor;
            uint64_t dist_prev;
        };

        struct par_gap_t {
            uint64_t beg;
            uint64_t end;
            uint64_t ref_beg;
            uint64_t dist_prev;
            uint64_t dist_next;
        };

        class merging_output {
        public:
            merging_output(factor_sink& sink, uint64_t& fact_counter)
                : output(sink)
                , num_fact(fact_counter)
            { }

            void add(uint64_t pos, factor f)
            {
                if (f.len > 0 && pending.len > 0 && pos == pending_pos + pending.len &&
                    pos - f.src == pending_pos - pending.src) {
                    pending.len = pending.len + f.len;
                    return;
                }

                flush();

                if (f.is_literal()) {
                    output(f);
                    num_fact++;
                    return;
                }

                pending = f;
                pending_pos = pos;
            }

            void flush()
            {
                if (pending.len == 0) return;
                output(pending);
                num_fact++;
                pending.len = 0;
            }

        private:
            factor_sink& output;
            uint64_t& num_fact;
            factor pending { .src = 0, .len = 0 };
            uint64_t pending_pos = 0;
        };

        struct par_buf_t {
            uint64_t num_pos = 0;
            std::vector<par_gap_t> gaps;
            std::vector<uint32_t> slot_refs;
            std::vector<uint32_t> slots;
            std::vector<uint40_t> vals;
            std::vector<uint64_t> owner_beg;
            std::vector<uint64_t> owner_fill;
        };

        factorizer(const text_t& input, parameters params)
            : time_start(params.time_start)
            , log(params.log)
            , p(params.num_threads)
            , tau(params.tau)
            , fact_mode(params.fact_mode)
            , transf_mode(params.transf_mode)
            , kind(params.range_ds)
            , exact_alg(params.exact_alg)
            , log_factor_count(params.log_factor_count)
            , T(input)
            , n(input.size())
        { }

        void factorize(quality_mode qual_mode, factor_sink& output)
        {
            if (n == 0) {
                return;
            }

            if (p == 0) {
                p = omp_get_max_threads();
            }

            baseline_bytes = malloc_count_current();
            malloc_count_reset_peak();

            struct omp_threads_guard {
                int saved = omp_get_max_threads();
                ~omp_threads_guard() { omp_set_num_threads(saved); }
            } omp_threads;

            omp_set_num_threads(p);

            if (log) {
                if (result_log::write_rows && result_log::out.is_open()) {
                    uint16_t transf_mode_int = qual_mode != exact ? 0 : (transf_mode + 1);

                    result_log::out << "RESULT"
                        << " text_name=" << result_log::text_name
                        << " n=" << n
                        << " alg=lz77_sss"
                        << " num_threads=" << p
                        << " tau=" << tau
                        << " fact_mode=" << fact_mode
                        << " transf_mode=" << transf_mode_int;
                }

                time = now();
                if (time_start.time_since_epoch().count() == 0) time_start = time;
            }

            if (qual_mode == exact && (exact_alg == sa_based && sa_supported())) {
                factorize_exact_sa(output);
            } else if (qual_mode == exact) {
                std::string aprx_file_name = std::filesystem::temp_directory_path().string()
                    + "/aprx_" + random_alphanumeric_string(10);
                direct_ofstream aprx_ofile(aprx_file_name);
                const uint8_t src_bytes = src_bytes_of(T);
                const uint8_t len_bytes = len_bytes_of(T);
                auto write_aprx = [&](factor f) { f.write(aprx_ofile, src_bytes, len_bytes); };
                factor_sink aprx_sink(write_aprx);
                compute_approximation(aprx_sink);
                aprx_sink.flush();
                aprx_ofile.close();
                const double peak_sa = estimated_sa_peak(n, sa_bytes(), T.size_in_bytes() + sa_extra_bytes());
                const double peak_sss = text_t::is_byte_text ? estimated_sss_peak(n, num_fact, transf_mode)
                    : estimated_sss_peak(n, T.size_in_bytes(), num_fact, transf_mode);

                if (exact_alg == auto_select && peak_sa <= max_sa_peak_ratio * peak_sss && sa_supported()) {
                    std::filesystem::remove(aprx_file_name);
                    LCE = lce_t();
                    release_free_memory();
                    factorize_exact_sa(output);
                } else {
                    uint64_t delta = std::min<uint64_t>(n / num_fact, max_delta);

                    exact_transformer(
                        T, n, LCE, aprx_file_name, delta, num_fact, p, log, transf_mode, kind)
                        .transform_to_exact(output);

                    std::filesystem::remove(aprx_file_name);
                }
            } else {
                compute_approximation(output);
            }

            if (log && fact_mode != skip_gaps) {
                uint64_t time_total = time_diff_ns(time_start, now());
                uint64_t peak_bytes = malloc_count_peak() - baseline_bytes + T.size_in_bytes();
                double comp_ratio = n / (double) num_fact;

                if (log_factor_count) {
                    std::cout << "num. of factors = " << num_fact << std::endl;
                    std::cout << "input length / num. of factors = " << comp_ratio << std::endl;
                }

                std::cout << "total time = " << format_time(time_total) << std::endl;
                std::cout << "throughput = " << format_throughput(n, time_total) << std::endl;
                std::cout << "peak memory consumption = " << format_size(peak_bytes);
                std::cout << " (" << (100.0 * peak_bytes) / n << " % of input)" << std::endl;

                if (result_log::write_rows && result_log::out.is_open()) {
                    result_log::out
                        << " num_factors=" << num_fact
                        << " comp_ratio=" << comp_ratio
                        << " time=" << time_total
                        << " throughput=" << throughput_mb_per_s(n, time_total)
                        << " mem_peak=" << peak_bytes << std::endl;
                }
            }
        }

        void compute_approximation(factor_sink& output)
        {

            build_LCE();
            build_LPF();
            LCE.free_sa_s();

            if (fact_mode == skip_gaps) {
                LPF.back().emplace_back(lpf_phrase { .beg = n, .end = n + 1, .src = 0 });
            } else {
                if (log) std::cout << "computing LPF statistics" << std::flush;

                compute_lpf_stats();
                LPF.back().emplace_back(lpf_phrase { .beg = n, .end = n + 1, .src = 0 });

                len_gaps = n - len_lpf_phr;
                double lpf_phr_per_sync = num_lpf_phr / (double)size_sss;
                double rel_len_gaps = len_gaps / (double)n;
                double gaps_per_lpf_phr = num_gaps / (double)num_lpf_phr;
                double avg_gap_len = len_gaps / (double)num_gaps;
                double avg_lpf_phr_len = len_lpf_phr / (double)num_lpf_phr;
                target_rh_idx_bytes = target_rh_idx_bytes_for(n, rel_len_gaps);
                double patt_len_guess = get_patt_len_guess(avg_gap_len, avg_lpf_phr_len, rel_len_gaps);
                const uint64_t cur_bytes = malloc_count_current() - baseline_bytes;
                const uint64_t budget_bytes = std::max<uint64_t>(malloc_count_peak() - baseline_bytes,
                    cur_bytes + rh_idx_t::size_in_bytes_for(n, target_rh_idx_bytes));
                auto exact_gaps_peak_bytes = [&](uint64_t ctx) {
                    return cur_bytes + exact_gaps_bytes(len_gaps + num_gaps * (1 + 2 * ctx), num_gaps, symbol_bytes);
                };
                factorize_gaps_exact = fact_mode == exact_gaps || exact_gaps_peak_bytes(0) <= budget_bytes;
                gap_ctx = 0;

                for (uint64_t step = max_gap_ctx_per_tau * tau; factorize_gaps_exact && step > 0; step /= 2) {
                    if (gap_ctx + step <= max_gap_ctx_per_tau * tau && exact_gaps_peak_bytes(gap_ctx + step) <= budget_bytes) {
                        gap_ctx += step;
                    }
                }

                if (!factorize_gaps_exact) {
                    gap_ctx = std::clamp<uint64_t>(hash_gap_ctx_factor / std::max(rel_len_gaps, 1e-9),
                        std::min<uint64_t>(min_hash_gap_ctx, tau), tau);
                }

                if (log) {
                    record_phase_time("lpf_stats", time_diff_ns(time, now()));
                    time = log_runtime(time);
                    std::cout << "|S| / (2n / tau) = " << size_sss / ((2.0 * n) / tau) << std::endl;
                    std::cout << "the density condition has " << (LCE.has_runs() ? "" : "not ") << "been applied" << std::endl;
                    std::cout << "num. of LPF phrases / SSS size = " << lpf_phr_per_sync << std::endl;
                    std::cout << "gaps length / input length = " << rel_len_gaps << std::endl;
                    std::cout << "num. of gaps / num. of LPF phrases = " << gaps_per_lpf_phr << std::endl;
                    std::cout << "avg. gap length = " << avg_gap_len << std::endl;
                    std::cout << "avg. LPF phrase length = " << avg_lpf_phr_len << std::endl;
                    std::cout << "pattern length guess = " << patt_len_guess << std::endl;
                    std::cout << "peak memory consumption = " << format_size(malloc_count_peak() - baseline_bytes + T.size_in_bytes()) << std::endl;
                    std::cout << "current memory consumption = " << format_size(malloc_count_current() - baseline_bytes + T.size_in_bytes()) << std::endl;
                    std::cout << "target index size = " << format_size(target_rh_idx_bytes) << std::endl;
                }

                if (!factorize_gaps_exact) {
                    for (auto [threshold, lens] : patt_len_table) {
                        if (patt_len_guess <= threshold) {
                            patt_lens = lens;
                            break;
                        }
                    }

                    if (log) {
                        std::cout << "pattern lengths for the rolling hash index: ";
                        for (uint64_t i = 0; i < num_patt_lens - 1; i++) std::cout << patt_lens[i] << ", ";
                        std::cout << patt_lens[num_patt_lens - 1] << std::endl;
                        std::cout << "initializing rolling hash index" << std::flush;
                    }

                    rh_idx = rh_idx_t(T, n, patt_lens, target_rh_idx_bytes, p);
                    if (log) std::cout << " (size = " << format_size(rh_idx.size_in_bytes()) << ")";

                    if (log) {
                        record_phase_time("init_rh_idx", time_diff_ns(time, now()));
                        time = log_runtime(time);
                    }
                }
            }

            factorize_from_lpf(output);
            LPF.clear();
            LPF.shrink_to_fit();
            rh_idx = rh_idx_t();
            release_free_memory();
        }

        inline uint64_t LCE_R(uint64_t i, uint64_t j) { return LCE.lce(i, j); }

        inline uint64_t LCE_L(uint64_t i, uint64_t j, uint64_t max_lce = std::numeric_limits<uint64_t>::max())
        {
            return T.lce_left(i, j, max_lce);
        }

        void build_LCE();

        void build_LPF();

        void compute_lpf_stats();

        inline factor longest_prev_occ(const par_buf_t& buf, uint64_t& gap, uint64_t pos);

        template <typename next_lpf_t>
        void prepare_block(next_lpf_t next_lpf, const par_blk_t& blk, uint64_t end,
            par_buf_t& buf, uint64_t num_owners, std::vector<uint64_t>& prefix_fps);

        void exchange_block(par_buf_t& buf);

        void update_rh_idx(std::vector<par_buf_t>& bufs, uint64_t owner, uint64_t num_blks_in_round);

        void output_blocks(merging_output& output, const std::vector<par_blk_t>& blks,
            const std::vector<std::array<std::vector<factor>, par_fact_bufs>>& facts,
            uint64_t round, uint64_t blks_per_round, uint64_t num_blks_in_round, uint64_t& out_end);

        template <typename next_lpf_t>
        void factorize_block(next_lpf_t next_lpf, const par_blk_t& blk, uint64_t end,
            const par_buf_t& buf, std::vector<factor>& facts);

        void factorize_from_lpf(factor_sink& output);

        void factorize_skip_gaps(factor_sink& output);

        template <typename lpf_beg_t, typename next_lpf_t>
        void factorize_hashed_gaps(factor_sink& output, lpf_beg_t lpf_beg, next_lpf_t next_lpf);

        template <typename lpf_beg_t, typename next_lpf_t>
        void factorize_exact_gaps(factor_sink& output, lpf_beg_t lpf_beg, next_lpf_t next_lpf);

        template <typename sym_t>
        static std::vector<sym_t> symbols_of(const text_t& text, uint16_t p);

        bool sa_supported() const;

        uint64_t sa_bytes() const;

        uint64_t sa_extra_bytes() const;

        void factorize_exact_sa(factor_sink& output);

        template <typename sa_t>
        void factorize_exact_sa(factor_sink& output);

        class exact_transformer {
        public:
            using sample_index_t = sample_index<text_t, lce_t, bit_aligned_vector>;
            using point_t = range_ds::point_t;
            using interval_t = sample_index_t::interval_t;
            using query_ctx_t = sample_index_t::query_ctx_t;

            std::chrono::steady_clock::time_point time_start, time;
            std::string aprx_file_name;
            std::string fact_file_name;
            bool log = false;
            bool write_queries = result_log::queries_out.is_open();
            uint16_t p = 0;
            transform_mode transf_mode = default_transf_mode;
            range_ds_kind kind;

            text_t T;
            const lce_t& LCE;

            uint64_t n = 0;
            uint64_t c = 0;
            uint64_t delta = 0;
            uint64_t num_aprx_fact = 0;
            uint64_t& num_fact;
            uint64_t num_par_sect;
            uint8_t src_bytes = 0;
            uint8_t len_bytes = 0;

            struct sect_info_t {
                uint64_t beg;
                uint64_t first_aprx;
            };

            std::vector<sect_info_t> par_sect;

            bit_aligned_vector C;
            sample_index_t idx_C;
            bit_aligned_interleaved_vectors<3> P;
            range_ds* R = nullptr;
            bit_aligned_vector Pi;
            bit_aligned_vector Psi;
            std::vector<uint8_t> is_smpl_len_left;

            inline uint64_t LCE_R(uint64_t i, uint64_t j) { return LCE.lce(i, j); }

            exact_transformer(const text_t& T, uint64_t n, const lce_t& LCE, std::string aprx_file_name,
                uint64_t delta, uint64_t& num_fact, uint16_t p, bool log, transform_mode transf_mode,
                range_ds_kind kind
            ) : aprx_file_name(aprx_file_name), log(log), p(p), transf_mode(transf_mode),
                kind(kind),
                T(T), LCE(LCE), n(n), delta(delta), num_aprx_fact(num_fact), num_fact(num_fact),
                src_bytes(src_bytes_of(T)), len_bytes(len_bytes_of(T)) { }

            exact_transformer(const exact_transformer&) = delete;
            exact_transformer& operator=(const exact_transformer&) = delete;

            ~exact_transformer() { delete R; }

            void transform_to_exact(factor_sink& output)
            {
                if (log) {
                    time = now();
                    time_start = time;
                }

                build_C();
                build_idx_C();
                build_P();
                build_Pi_Psi();

                if (log) {
                    std::cout << "building " << kind.name() << std::flush;
                }

                if (kind.decomposed) {
                    R = make_range_ds(kind, T, C, P, p);

                    if (kind.is_static() && !write_queries) P.clear();
                } else {
                    std::vector<point_t> points;
                    no_init_resize(points, c);

                    #pragma omp parallel for num_threads(p) schedule(dynamic, 65536)
                    for (uint64_t i = 0; i < c; i++) {
                        points[i] = point_t { .x = P.get<0>(i), .y = P.get<1>(i),
                            .weight = P.get<2>(i) };
                    }

                    if (kind.is_static() && !write_queries) P.clear();

                    R = make_range_ds(kind, points, c, p);
                }

                if (log) {
                    record_phase_time("range_ds", time_diff_ns(time, now()));
                    std::cout << " (" << format_size(R->size_in_bytes()) << ")";
                    time = log_runtime(time);
                }

                if (R->is_dynamic()) {
                    num_par_sect = 1;
                    par_sect.resize(2);
                    par_sect[1] = {.beg = n, .first_aprx = num_aprx_fact};
                }

                if (p > 1) {
                    fact_file_name = std::filesystem::temp_directory_path().string()
                        + "/fact_" + random_alphanumeric_string(10);
                }

                if (transf_mode == with_interval_samples) {
                    transform_to_exact_with_interval_samples(output);
                } else if (transf_mode == without_interval_samples) {
                    transform_to_exact_without_interval_samples(output);
                }

                if (log) {
                    record_phase_time("compute_exact", time_diff_ns(time, now()));
                    time = log_runtime(time);
                }

                if (log && R->is_dynamic()) {
                    std::cout << "final size of " << R->name()
                              << " = " << format_size(R->size_in_bytes()) << std::endl;
                }
            }

            void build_C();

            void build_idx_C();

            void build_Pi_Psi();

            void build_P();

            void insert_points_before(uint64_t& x_r, uint64_t i);

            void find_close_sources(factor& f, uint64_t i, uint64_t e);

            inline void seek_x_c(uint64_t& x_c, uint64_t pos);

            bool intersect(
                const interval_t& pa_c_iv, const interval_t& sa_c_iv,
                uint64_t i, uint64_t j, uint64_t lce_l, uint64_t lce_r, uint64_t& x_c, factor& f);

            void improve_factor_without_interval_samples(uint64_t i, uint64_t e, uint64_t& x_c, factor& f);

            void transform_to_exact_without_interval_samples(factor_sink& output);

            void extend_right_with_interval_samples(
                const interval_t& pa_c_iv,
                uint64_t i, uint64_t j, uint64_t e, uint64_t& x_c, factor& f);

            void improve_factor_with_interval_samples(
                uint64_t i, uint64_t e, uint64_t& x_c, std::vector<uint64_t>& fp_left, factor& f);

            void transform_to_exact_with_interval_samples(factor_sink& output);

            factor aprx_factor_at(uint64_t i);

            template <typename factor_at_t>
            void combine_factorizations(factor_sink& output, factor_at_t factor_at);
        };
    };
};

#ifndef LZ77_SSS_DECLARATIONS_ONLY
#include "algorithms/common.cpp"

#include "algorithms/aprx/lpf_stats.cpp"
#include "algorithms/aprx/fact/common.cpp"
#include "algorithms/aprx/fact/hashed_gaps.cpp"
#include "algorithms/aprx/fact/exact_gaps.cpp"
#include "algorithms/aprx/fact/skip_gaps.cpp"
#include "algorithms/aprx/lpf.cpp"

#include "algorithms/aprx_to_exact/common.cpp"
#include "algorithms/aprx_to_exact/with_interval_samples.cpp"
#include "algorithms/aprx_to_exact/without_interval_samples.cpp"

#include "algorithms/exact_sa.cpp"
#endif