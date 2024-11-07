#pragma once

#include <filesystem>
#include <functional>
#include <vector>

#include <lce/lce_classic.hpp>
#include <lce/lce_naive_wordwise_xor.hpp>
#include <lce/lce_sss.hpp>

#include <lz77_sss/algorithms/lce_l.hpp>
#include <lz77_sss/data_structures/dynamic_range/dynamic_square_grid.hpp>
#include <lz77_sss/data_structures/dynamic_range/semi_dynamic_square_grid.hpp>
#include <lz77_sss/data_structures/parallel_rolling_hash_index_107.hpp>
#include <lz77_sss/data_structures/rolling_hash_index_107.hpp>
#include <lz77_sss/data_structures/rolling_hash_index_61.hpp>
#include <lz77_sss/data_structures/sample_index/sample_index.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_kd_tree.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_square_grid.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_striped_square.hpp>
#include <lz77_sss/misc/utils.hpp>

enum phrase_mode {
    lpf_naive,
    lpf_lnf_naive,
    lpf_opt,
    lpf_lnf_opt
};

enum factorize_mode {
    greedy_naive,
    greedy,
    skip_phrases
};

enum transform_mode {
    naive,
    with_samples,
    without_samples,
};

struct parameters {
    uint16_t num_threads = 0;
    bool log = false;
};

template <typename pos_t = uint32_t>
class lz77_sss {
public:
    static_assert(std::is_same_v<pos_t, uint32_t> || std::is_same_v<pos_t, uint64_t>);

    static constexpr phrase_mode        default_phr_mode         = lpf_opt;
    static constexpr factorize_mode     default_fact_mode        = greedy;
    static constexpr transform_mode     default_transf_mode      = without_samples;
    template <typename sidx_t> using    default_range_ds_t       = decomposed_static_weighted_kd_tree<sidx_t>;

    static constexpr uint64_t           default_tau              = 512;
    static constexpr uint64_t           max_delta                = 256;
    static constexpr uint64_t           range_scan_threshold     = 4096;
    static constexpr uint64_t           max_buffered_lpf_phrases = 2048;
    static constexpr double             min_rel_gap_len          = 0.2;
    static constexpr uint64_t           min_gap_blk_size         = 4096;
    static constexpr uint64_t           max_num_gap_blks         = 512;
    static constexpr uint64_t           num_patt_lens            = 5;

    using entry_t = std::pair<pos_t, std::array<pos_t, num_patt_lens>>;

    static constexpr std::array<entry_t, 10> patt_len_table {
        entry_t { 6, { 2, 3, 4, 5, 6 } },
        entry_t { 8, { 2, 3, 4, 6, 8 } },
        entry_t { 12, { 2, 3, 4, 8, 12 } },
        entry_t { 16, { 2, 4, 6, 9, 16 } },
        entry_t { 32, { 2, 4, 6, 10, 20 } },
        entry_t { 64, { 2, 4, 7, 12, 28 } },
        entry_t { 128, { 2, 4, 8, 16, 36 } },
        entry_t { 256, { 2, 5, 10, 20, 42 } },
        entry_t { 1024, { 2, 6, 12, 24, 48 } },
        entry_t { UINT_MAX, { 2, 8, 16, 32, 64 } }
    };
    
    static double get_patt_len_guess(double avg_gap_len, double avg_lpf_phr_len, double rel_len_gaps)
    {
        return std::min<double>({avg_gap_len, avg_lpf_phr_len, 8.0 * std::pow(128, 1.0 - rel_len_gaps)});
    }

    static uint64_t get_target_gap_idx_size(uint64_t n, double rel_len_gaps)
    {
        return std::max<uint64_t>(malloc_count_peak() - malloc_count_current(), (n / 3.0) * rel_len_gaps);
    }

    struct factor {
        pos_t src;
        pos_t len;

        friend class lz77_sss;

        friend std::istream& operator>>(std::istream& in, const factor& f)
        {
            in.read((char*) &f, sizeof(factor));
            return in;
        }

        friend std::ostream& operator<<(std::ostream& out, const factor& f)
        {
            out.write((char*) &f, sizeof(factor));
            return out;
        }
    };

    using output_it_t = std::function<void(factor)>;

    template <
        factorize_mode                  fact_mode   = default_fact_mode,
        phrase_mode                     phr_mode    = default_phr_mode,
        uint64_t                        tau         = default_tau
    >
    static void factorize_approximate(std::string& input, output_it_t output, parameters params = { })
    {
        factorizer<tau>(input, params).template factorize<approximate, fact_mode, phr_mode>(output);
    }

    template <
        factorize_mode                  fact_mode   = default_fact_mode,
        phrase_mode                     phr_mode    = default_phr_mode,
        uint64_t                        tau         = default_tau,
        std::output_iterator<factor>    out_it_t
    >
    static void factorize_approximate(std::string& input, out_it_t output, parameters params = { })
    {
        factorize_approximate<fact_mode, phr_mode, tau>(input, [&](factor f){*output++ = f;}, params);
    }

    template <
        factorize_mode  fact_mode   = default_fact_mode,
        phrase_mode     phr_mode    = default_phr_mode,
        uint64_t        tau         = default_tau
    >
    static void factorize_approximate(std::string& input, std::ofstream& out, parameters params = { })
    {
        factorize_approximate<fact_mode, phr_mode, tau>(input, std::ostream_iterator<factor>(out, ""), params);
    }

    template <
        factorize_mode  fact_mode   = default_fact_mode,
        phrase_mode     phr_mode    = default_phr_mode,
        uint64_t        tau         = default_tau
    >
    static std::vector<factor> factorize_approximate(std::string& input, parameters params = { })
    {
        std::vector<factor> factorization;
        factorize_approximate<fact_mode, phr_mode, tau>(input, std::back_insert_iterator(factorization), params);
        return factorization;
    }

    template <
        factorize_mode                  fact_mode       = default_fact_mode,
        phrase_mode                     phr_mode        = default_phr_mode,
        transform_mode                  transf_mode     = default_transf_mode,
        template <typename> typename    range_ds_t      = default_range_ds_t,
        uint64_t                        tau             = default_tau
    >
    static void factorize_exact(std::string& input, output_it_t output, parameters params = { })
    {
        factorizer<tau>(input, params).template factorize<exact, fact_mode, phr_mode, transf_mode, range_ds_t>(output);
    }

    template <
        factorize_mode                  fact_mode       = default_fact_mode,
        phrase_mode                     phr_mode        = default_phr_mode,
        transform_mode                  transf_mode     = default_transf_mode,
        template <typename> typename    range_ds_t      = default_range_ds_t,
        uint64_t                        tau             = default_tau,
        std::output_iterator<factor>    out_it_t
    >
    static void factorize_exact(std::string& input, out_it_t output, parameters params = { })
    {
        factorize_exact<fact_mode, phr_mode, transf_mode, range_ds_t, tau>(input, [&](factor f){*output++ = f;}, params);
    }

    template <
        factorize_mode                  fact_mode       = default_fact_mode,
        phrase_mode                     phr_mode        = default_phr_mode,
        transform_mode                  transf_mode     = default_transf_mode,
        template <typename> typename    range_ds_t      = default_range_ds_t,
        uint64_t                        tau             = default_tau
    >
    static void factorize_exact(std::string& input, std::ofstream& out, parameters params = { })
    {
        factorize_exact<fact_mode, phr_mode, transf_mode, range_ds_t, tau>(
            input, std::ostream_iterator<factor>(out, ""), params);
    }

    template <
        factorize_mode                  fact_mode       = default_fact_mode,
        phrase_mode                     phr_mode        = default_phr_mode,
        transform_mode                  transf_mode     = default_transf_mode,
        template <typename> typename    range_ds_t      = default_range_ds_t,
        uint64_t                        tau             = default_tau
    >
    static std::vector<factor> factorize_exact(std::string& input, parameters params = { })
    {
        std::vector<factor> factorization;
        factorize_exact<fact_mode, phr_mode, transf_mode, range_ds_t, tau>(
            input, std::back_insert_iterator(factorization), params);
        return factorization;
    }
    
    template <std::input_iterator fact_it_t>
    static std::string decode(fact_it_t fact_it, pos_t output_size);

    static std::string decode(std::vector<factor>& factorization, pos_t output_size)
    {
        return decode(factorization.begin(), output_size);
    }

    static std::string decode(std::ifstream& in, pos_t output_size)
    {
        return decode(std::istream_iterator<factor>(in), output_size);
    }

protected:
    struct lpf {
        pos_t beg;
        pos_t end;
        pos_t src;

        friend std::istream& operator>>(std::istream& in, const lpf& p)
        {
            in.read((char*) &p, sizeof(lpf));
            return in;
        }

        friend std::ostream& operator<<(std::ostream& out, const lpf& p)
        {
            out.write((char*) &p, sizeof(lpf));
            return out;
        }
    };

    enum quality_mode {
        approximate,
        exact
    };

    lz77_sss() = delete;
    lz77_sss(lz77_sss&& other) = delete;
    lz77_sss(const lz77_sss& other) = delete;
    lz77_sss& operator=(lz77_sss&& other) = delete;
    lz77_sss& operator=(const lz77_sss& other) = delete;

    template <uint64_t tau>
    class factorizer {
    public:
        using lce_t = alx::lce::lce_sss<char, tau, pos_t, false>;
        using gap_idx_t = rolling_hash_index_107<pos_t, num_patt_lens>;
        using par_gap_idx_t = parallel_rolling_hash_index_107<pos_t, num_patt_lens>;
        using fp_arr_t = par_gap_idx_t::fp_arr_t;
        std::chrono::steady_clock::time_point time_start, time, time_end;
        uint64_t baseline_memory_alloc = 0;
        uint64_t target_index_size = 0;
        bool log = false;
        uint16_t p = 0;
        pos_t roll_threshold = 0;

        std::string& T;
        pos_t n = 0;
        pos_t size_sss = 0;
        pos_t num_lpf = 0;
        pos_t len_lpf_phr = 0;
        pos_t num_phr = 0;
        pos_t len_gaps = 0;
        pos_t num_gaps = 0;

        lce_t LCE;
        std::vector<std::vector<lpf>> LPF;
        std::array<pos_t, num_patt_lens> patt_lens;
        gap_idx_t gap_idx;
        par_gap_idx_t par_gap_idx;

        std::vector<uint32_t> PSV_S;
        std::vector<uint32_t> NSV_S;
        std::vector<uint32_t> PGV_S;
        std::vector<uint32_t> NGV_S;

        struct lpf_arr_it_t {
            uint16_t i_p;
            uint32_t i;
        };
        
        struct block_info_t {
            uint16_t i_p;
            uint32_t i;
            pos_t beg;
        };

        std::vector<block_info_t> blk_info;
        std::vector<std::vector<factor>> factors;

        factorizer(std::string& input, parameters params)
            : log(params.log)
            , p(params.num_threads)
            , T(input)
            , n(input.size())
        { }

        template <
            quality_mode qual_mode,
            factorize_mode fact_mode = default_fact_mode,
            phrase_mode phr_mode = default_phr_mode,
            transform_mode transf_mode = default_transf_mode,
            template <typename> typename range_ds_t = default_range_ds_t>
        void factorize(output_it_t& output)
        {
            if (log) {
                #ifdef LZ77_SSS_BENCH
                if (result_file_path != "") {
                    uint16_t transf_mode_int = qual_mode != exact ? 0 : (transf_mode + 1);

                    result_file << "RESULT"
                        << " text_name=" << text_name
                        << " n=" << n
                        << " alg=lz77_sss"
                        << " num_threads=" << p
                        << " tau=" << tau
                        << " phr_mode=" << phr_mode
                        << " fact_mode=" << fact_mode
                        << " transf_mode=" << transf_mode_int;
                }
                #endif

                time = now();
                time_start = time;
            }

            if (p == 0) {
                p = omp_get_max_threads();
            }

            baseline_memory_alloc = malloc_count_current();
            malloc_count_reset_peak();
            omp_set_num_threads(p);

            if constexpr (qual_mode == exact) {
                static_assert(fact_mode != skip_phrases);
                std::string aprx_file_name = std::filesystem::temp_directory_path().string()
                    + "/aprx_" + random_alphanumeric_string(10);
                std::ofstream ofile_aprx(aprx_file_name);
                std::ostream_iterator<factor> ofile_aprx_it(ofile_aprx, "");
                output_it_t aprx_out_it = [&](factor f) { *ofile_aprx_it++ = f; };
                compute_approximation<fact_mode, phr_mode>(aprx_out_it);
                ofile_aprx.close();
                std::ifstream ifile_aprx(aprx_file_name);
                std::istream_iterator<factor> ifile_aprx_it(ifile_aprx);
                pos_t delta = std::min<pos_t>(n / num_phr, max_delta);

                if constexpr (std::is_same_v<pos_t, uint32_t>) {
                    exact_factorizer<uint32_t, transf_mode, range_ds_t>(
                        T, LCE, delta, num_phr, p, log)
                        .transform_to_exact(ifile_aprx_it, output);
                } else {
                    pos_t max_num_samples = num_phr + n / delta;

                    if (max_num_samples <= std::numeric_limits<uint32_t>::max()) {
                        exact_factorizer<uint32_t, transf_mode, range_ds_t>(
                            T, LCE, delta, num_phr, p, log)
                            .transform_to_exact(ifile_aprx_it, output);
                    } else {
                        exact_factorizer<uint64_t, transf_mode, range_ds_t>(
                            T, LCE, delta, num_phr, p, log)
                            .transform_to_exact(ifile_aprx_it, output);
                    }
                }

                ifile_aprx.close();
                std::filesystem::remove(aprx_file_name);
            } else {
                compute_approximation<fact_mode, phr_mode>(output);
            }

            if (log && fact_mode != skip_phrases) {
                uint64_t time_total = time_diff_ns(time_start, now());
                uint64_t mem_peak = malloc_count_peak() - baseline_memory_alloc;
                double comp_ratio = n / (double)num_phr;

                std::cout << "compression ratio: " << comp_ratio << std::endl;
                std::cout << "total time: " << format_time(time_total) << std::endl;
                std::cout << "throughput: " << format_throughput(n, time_total) << std::endl;
                std::cout << "peak memory consumption: " << format_size(mem_peak) << std::endl;

                #ifdef LZ77_SSS_BENCH
                if (result_file_path != "") {
                    result_file
                        << " num_factors=" << num_phr
                        << " comp_ratio=" << comp_ratio
                        << " time=" << time_total
                        << " throughput=" << throughput(n, time_total)
                        << " mem_peak=" << mem_peak << std::endl;
                }
                #endif
            }
        }

        template <
            factorize_mode fact_mode = default_fact_mode,
            phrase_mode phr_mode = default_phr_mode>
        void compute_approximation(output_it_t& output)
        {
            LPF.resize(p);

            if constexpr (phr_mode == lpf_naive) {
                build_lce();
                build_LPF_naive();
            } else if constexpr (phr_mode == lpf_opt) {
                build_lce();
                build_LPF_opt([&](uint16_t i_p, lpf&& p) {
                    LPF[i_p].emplace_back(p);});
            } else {
                if (log) std::cout << "reversing T" << std::flush;
                std::reverse(T.begin(), T.end());
                if (log) time = log_runtime(time);
                build_lce();
                build_LNF_all<phr_mode>([&](uint16_t i_p, lpf&& p) {
                    LPF[i_p].emplace_back(p);});
                LCE = lce_t();
                if (log) std::cout << "reversing T" << std::flush;
                std::reverse(T.begin(), T.end());
                if (log) time = log_runtime(time);
                build_lce();
                build_LPF_all<phr_mode>([&](uint16_t i_p, lpf&& p) {
                    LPF[i_p].emplace_back(p);});
            }

            LCE.delete_ssa();

            if constexpr (phr_mode == lpf_lnf_naive || phr_mode == lpf_lnf_opt) {
                if (log) {
                    std::cout << "selecting LPF phrases" << std::flush;
                }

                #pragma omp parallel num_threads(p)
                {
                    uint16_t i_p = omp_get_thread_num();
                    greedy_phrase_selection(LPF[i_p]);
                }

                if (log) {
                    log_phase("select_phrases", time_diff_ns(time, now()));
                    time = log_runtime(time);
                }
            }

            if constexpr (fact_mode == skip_phrases) {
                LPF[p - 1].emplace_back(lpf { .beg = n, .end = n + 1 });
            } else {
                if (log) std::cout << "computing LPF statistics" << std::flush;

                get_phrase_info();
                LPF[p - 1].emplace_back(lpf { .beg = n, .end = n + 1 });

                len_gaps = n - len_lpf_phr;
                double lpf_phr_per_sync = num_lpf / (double)size_sss;
                double rel_len_gaps = len_gaps / (double)n;
                double gaps_per_lpf_phr = num_gaps / (double)num_lpf;
                double avg_gap_len = len_gaps / (double)num_gaps;
                double avg_lpf_phr_len = len_lpf_phr / (double)num_lpf;
                target_index_size = get_target_gap_idx_size(n, rel_len_gaps);

                if (log) {
                    log_phase("phrase_info", time_diff_ns(time, now()));
                    time = log_runtime(time);
                    std::cout << "num. of LPF phrases / SSS size = " << lpf_phr_per_sync << std::endl;
                    std::cout << "gaps length / input length: " << rel_len_gaps << std::endl;
                    std::cout << "num. of gaps / num. of LPF phrases: " << gaps_per_lpf_phr << std::endl;
                    std::cout << "avg. gap length: " << avg_gap_len << std::endl;
                    std::cout << "avg. LPF phrase length: " << avg_lpf_phr_len << std::endl;
                    std::cout << "peak memory consumption: " << format_size(malloc_count_peak() - baseline_memory_alloc) << std::endl;
                    std::cout << "current memory consumption: " << format_size(malloc_count_current() - baseline_memory_alloc) << std::endl;
                    std::cout << "target index size: " << format_size(target_index_size) << std::endl;
                }

                double patt_len_guess = get_patt_len_guess(avg_gap_len, avg_lpf_phr_len, rel_len_gaps);

                for (auto [threshold, lens] : patt_len_table) {
                    if (patt_len_guess <= threshold) {
                        patt_lens = lens;
                        break;
                    }
                }

                for_constexpr<0, num_patt_lens, 1>([&](auto j) {
                    roll_threshold += patt_lens[j];
                });

                roll_threshold /= num_patt_lens;

                if (log) {
                    std::cout << "pattern lengths for the rolling hash index: ";
                    for (pos_t i = 0; i < num_patt_lens - 1; i++) std::cout << patt_lens[i] << ", ";
                    std::cout << patt_lens[num_patt_lens - 1] << std::endl;
                    std::cout << "initializing rolling hash index" << std::flush;
                }

                if (fact_mode == greedy && rel_len_gaps > min_rel_gap_len && p >= 2) {
                    par_gap_idx = par_gap_idx_t(T.data(), n, patt_lens, target_index_size, p);
                    if (log) std::cout << " (size: " << format_size(par_gap_idx.size_in_bytes()) << ")";
                } else {
                    gap_idx = gap_idx_t(T.data(), n, patt_lens, target_index_size);
                    if (log) std::cout << " (size: " << format_size(gap_idx.size_in_bytes()) << ")";
                }

                if (log) {
                    log_phase("init_gap_idx", time_diff_ns(time, now()));
                    time = log_runtime(time);
                }
            }

            factorize<fact_mode>(output);
            LPF.clear();
            LPF.shrink_to_fit();
            gap_idx = gap_idx_t();
            par_gap_idx = par_gap_idx_t();
        }

        inline pos_t LCE_R(pos_t i, pos_t j)
        {
            return LCE.lce(i, j);
        }

        inline pos_t LCE_L(pos_t i, pos_t j, pos_t max_lce = std::numeric_limits<pos_t>::max())
        {
            return lce_l_128<pos_t>(T.data(), i, j, max_lce);
        }

        static void greedy_phrase_selection(std::vector<lpf>& P);

        void build_lce();

        void build_PSV_NSV_S();

        void build_PGV_NGV_S();

        void build_LPF_naive();

        void build_LPF_opt(std::function<void(uint16_t, lpf&&)> lpf_it);

        template <phrase_mode phr_mode>
        void build_LNF_all(std::function<void(uint16_t, lpf&&)> lpf_it);
        
        template <phrase_mode phr_mode>
        void build_LPF_all(std::function<void(uint16_t, lpf&&)> lpf_it);

        void get_phrase_info();

        inline factor longest_prev_occ(pos_t pos);

        template <bool first_block>
        inline factor longest_prev_occ_par(fp_arr_t& fps, pos_t pos, pos_t blk_end);

        template <bool first_block, typename lpf_it_t>
        void factorize_block(std::function<lpf(lpf_it_t&)>& next_lpf, pos_t blk_beg, pos_t blk_end);

        template <factorize_mode fact_mode>
        void factorize(output_it_t& output);

        template <factorize_mode fact_mode, typename lpf_it_t>
        void factorize_sequential(output_it_t& output, std::function<lpf_it_t()>& lpf_beg, std::function<lpf(lpf_it_t&)>& next_lpf)
        {
            if constexpr (fact_mode == skip_phrases) {
                factorize_skip_gaps(output, lpf_beg, next_lpf);
            } else if constexpr (fact_mode == greedy) {
                factorize_greedy(output, lpf_beg, next_lpf);
            } else if constexpr (fact_mode == greedy_naive) {
                factorize_greedy_naive(output, lpf_beg, next_lpf);
            }
        }

        template <typename lpf_it_t>
        void factorize_skip_gaps(output_it_t& output, std::function<lpf_it_t()>& lpf_beg, std::function<lpf(lpf_it_t&)>& next_lpf);

        template <typename lpf_it_t>
        void factorize_greedy_naive(output_it_t& output, std::function<lpf_it_t()>& lpf_beg, std::function<lpf(lpf_it_t&)>& next_lpf);

        template <typename lpf_it_t>
        void factorize_greedy(output_it_t& output, std::function<lpf_it_t()>& lpf_beg, std::function<lpf(lpf_it_t&)>& next_lpf);

        template <typename lpf_it_t>
        void factorize_greedy_parallel(output_it_t& output, std::function<lpf_it_t()>& lpf_beg, std::function<lpf(lpf_it_t&)>& next_lpf);

        template <
            typename sidx_t,
            transform_mode transf_mode,
            template <typename> typename range_ds_t>
        class exact_factorizer {
        public:
            using sample_index_t = sample_index<pos_t, sidx_t, lce_t>;
            using point_t = typename range_ds_t<sidx_t>::point_t;
            using sxa_interval_t = sample_index<pos_t, sidx_t, lce_t>::sxa_interval_t;
            using query_context_t = sample_index<pos_t, sidx_t, lce_t>::query_context;

            std::chrono::steady_clock::time_point time_start, time, time_end;
            std::string fact_file_name;
            bool log = false;
            uint16_t p = 0;

            const std::string& T;
            const lce_t& LCE;

            pos_t n = 0;
            pos_t c = 0;
            pos_t delta = 0;
            pos_t& num_phr;

            std::vector<pos_t> start_thr;

            std::vector<pos_t> C;
            sample_index_t idx_C;
            std::vector<point_t> P;
            range_ds_t<sidx_t> R;
            std::vector<sidx_t> PS;
            std::vector<sidx_t> SP;

            inline pos_t LCE_R(pos_t i, pos_t j)
            {
                return LCE.lce(i, j);
            }

            exact_factorizer(std::string& T, const lce_t& LCE, pos_t delta, pos_t& num_phr, uint16_t p, bool log)
                : log(log), p(p), T(T), LCE(LCE), n(T.size()), delta(delta), num_phr(num_phr) { }

            void transform_to_exact(std::istream_iterator<factor>& ifile_aprx_it, output_it_t& output)
            {
                if (log) {
                    time = now();
                    time_start = time;
                }

                build_c(ifile_aprx_it);
                build_idx_C();
                build_p();

                if constexpr (transf_mode != naive) {
                    build_ps_sp();
                }

                if (log) {
                    std::cout << "building " << range_ds_t<sidx_t>::name() << std::flush;
                }

                if constexpr (range_ds_t<sidx_t>::is_decomposed()) {
                    R = range_ds_t<sidx_t>(T, C, P, p);
                } else {
                    R = range_ds_t<sidx_t>(P, c, p);
                }

                if constexpr (range_ds_t<sidx_t>::is_static()) {
                    #if not defined(GEN_RANGE_QUERIES)
                    P.clear();
                    P.shrink_to_fit();
                    #endif
                }

                if (log) {
                    log_phase("range_ds", time_diff_ns(time, now()));
                    std::cout << " (" << format_size(R.size_in_bytes()) << ")";
                    time = log_runtime(time);
                }

                if (range_ds_t<sidx_t>::is_dynamic()) {
                    p = 1;
                    start_thr.resize(p + 1);
                    start_thr[p] = n;
                }

                if (p > 1) {
                    fact_file_name = std::filesystem::temp_directory_path().string()
                        + "/fact_" + random_alphanumeric_string(10);
                }

                if constexpr (transf_mode == naive) {
                    transform_to_exact_naive(output);
                } else if constexpr (transf_mode == with_samples) {
                    transform_to_exact_with_samples(output);
                } else if constexpr (transf_mode == without_samples) {
                    transform_to_exact_without_samples(output);
                }

                if (log) {
                    log_phase("compute_exact", time_diff_ns(time, now()));
                    time = log_runtime(time);
                }

                if (log && range_ds_t<sidx_t>::is_dynamic()) {
                    std::cout << "final size of " << range_ds_t<sidx_t>::name()
                              << ": " << format_size(R.size_in_bytes()) << std::endl;
                }
            }

            void build_c(std::istream_iterator<factor>& ifile_aprx_it);

            void build_idx_C();

            void build_ps_sp();

            void build_p();

            void insert_points(sidx_t& x_c, pos_t i);

            void find_close_sources(factor& f, pos_t i, pos_t e);

            inline void adjust_xc(sidx_t& gap_idx, pos_t pos);

            bool intersect(
                const sxa_interval_t& spa_iv, const sxa_interval_t& ssa_iv,
                pos_t i, pos_t j, pos_t lce_l, pos_t lce_r, sidx_t& x_c, factor& f);

            void transform_to_exact_naive(output_it_t& output);

            void transform_to_exact_without_samples(output_it_t& output);

            void extend_right_with_samples(
                const sxa_interval_t& spa_iv,
                pos_t i, pos_t j, pos_t e, sidx_t& x_c, factor& f);

            void transform_to_exact_with_samples(output_it_t& output);

            void combine_factorizations(output_it_t& output);
        };
    };
};

#include "algorithms/common.cpp"

#include "algorithms/approximate/common.cpp"
#include "algorithms/approximate/factorize/common.cpp"
#include "algorithms/approximate/factorize/greedy.cpp"
#include "algorithms/approximate/factorize/greedy_naive.cpp"
#include "algorithms/approximate/factorize/greedy_parallel.cpp"
#include "algorithms/approximate/factorize/skip_gaps.cpp"
#include "algorithms/approximate/lpf_lnf/nxv_pxv.cpp"
#include "algorithms/approximate/lpf_lnf/lpf_naive.cpp"
#include "algorithms/approximate/lpf_lnf/lpf_opt.cpp"
#include "algorithms/approximate/lpf_lnf/lpf_lnf.cpp"

#include "algorithms/transform_to_exact/common.cpp"
#include "algorithms/transform_to_exact/naive.cpp"
#include "algorithms/transform_to_exact/with_samples.cpp"
#include "algorithms/transform_to_exact/without_samples.cpp"