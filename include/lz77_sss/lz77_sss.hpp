#pragma once

#include <filesystem>
#include <functional>
#include <vector>

#include <ds/lce_classic.hpp>
#include <ds/lce_naive_wordwise_xor.hpp>
#include <ds/lce_sss.hpp>

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
    template <typename sidx_t> using    default_range_ds_t       = decomposed_static_weighted_square_grid<sidx_t>;

    static constexpr uint64_t           default_tau              = 512;
    static constexpr uint64_t           max_delta                = 256;
    static constexpr uint64_t           rks_sample_rate          = 16;
    static constexpr uint64_t           range_scan_threshold     = 4096;
    static constexpr uint64_t           min_par_input_size       = 500'000;
    static constexpr double             min_par_rel_gap_len      = 0.2;
    static constexpr uint64_t           min_par_gap_blk_size     = 4096;
    static constexpr uint64_t           max_par_gap_blks         = 512;
    static constexpr uint64_t           num_patt_lens            = 5;
    static constexpr uint64_t           min_rh_index_size        = 1 << 20;
    static constexpr uint64_t           max_rh_index_size        = 1 << 30;
    static constexpr double             min_rel_rh_index_size    = 0.1;

    using entry_t = std::pair<double, std::array<pos_t, num_patt_lens>>;
    static constexpr double infty = std::numeric_limits<double>::max();

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
        entry_t { infty, { 2, 8, 16, 32, 64 } }
    };
    
    static double get_patt_len_guess(double avg_gap_len, double avg_lpf_phr_len, double rel_len_gaps)
    {
        return std::min<double>({avg_gap_len, avg_lpf_phr_len, 8.0 * std::pow(128, 1.0 - rel_len_gaps)});
    }

    static uint64_t get_target_gap_idx_size(uint64_t n, double rel_len_gaps)
    {
        return std::min<uint64_t>(max_rh_index_size, std::max<uint64_t>({min_rh_index_size,
            malloc_count_peak() - malloc_count_current(), (n / 3.0) * rel_len_gaps}));
    }

    static uint64_t get_max_smpl_len_right(double aprx_comp_ratio)
    {
        return std::round(aprx_comp_ratio * (1.0 + 0.5 * std::exp(-aprx_comp_ratio / 1000.0)));
    }

    struct factor {
        pos_t src;
        pos_t len;

        friend class lz77_sss;

        pos_t length() const
        {
            return std::max<pos_t>(1, len);
        }

        static constexpr pos_t size_of()
        {
            if constexpr (std::is_same_v<pos_t, uint32_t>) {
                return 8;
            } else {
                return 10;
            }
        }

        friend std::istream& operator>>(std::istream& in, factor& f)
        {
            if constexpr (std::is_same_v<pos_t, uint32_t>) {
                in.read((char*) &f, 8);
            } else {
                f.src = 0;
                f.len = 0;
                in.read((char*) &f.src, 5);
                in.read((char*) &f.len, 5);
            }

            return in;
        }

        friend std::ostream& operator<<(std::ostream& out, const factor& f)
        {
            if constexpr (std::is_same_v<pos_t, uint32_t>) {
                out.write((char*) &f, 8);
            } else {
                out.write((char*) &f.src, 5);
                out.write((char*) &f.len, 5);
            }

            return out;
        }
    };

    template <
        factorize_mode  fact_mode = default_fact_mode,
        phrase_mode     phr_mode  = default_phr_mode,
        uint64_t        tau       = default_tau,
        typename char_t,
        typename output_fnc_t
    >
    static void factorize_approximate(char_t* input, pos_t input_size, output_fnc_t output, parameters params = { })
    {
        factorizer<tau, char_t>(input, input_size, params).template factorize<approximate, fact_mode, phr_mode>(output);
    }

    template <
        factorize_mode               fact_mode   = default_fact_mode,
        phrase_mode                  phr_mode    = default_phr_mode,
        transform_mode               transf_mode = default_transf_mode,
        template <typename> typename range_ds_t  = default_range_ds_t,
        uint64_t                     tau         = default_tau,
        typename char_t,
        typename output_fnc_t
    >
    static void factorize_exact(char_t* input, pos_t input_size, output_fnc_t output, parameters params = { })
    {
        factorizer<tau, char_t>(input, input_size, params).template factorize<exact, fact_mode, phr_mode, transf_mode, range_ds_t>(output);
    }
    
    template <std::input_iterator fact_it_t, typename char_t>
    static void decode(fact_it_t fact_it, char_t* output, pos_t output_size);

protected:
    struct lpf {
        pos_t beg;
        pos_t end;
        pos_t src;
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

    template <uint64_t tau, typename char_t>
    class factorizer {
    public:
        using lce_t = lce::ds::lce_sss<char_t, tau, pos_t, false>;
        using gap_idx_t = rolling_hash_index_107<pos_t, num_patt_lens, char_t>;
        using par_gap_idx_t = parallel_rolling_hash_index_107<pos_t, num_patt_lens, char_t>;
        using fp_arr_t = par_gap_idx_t::fp_arr_t;
        std::chrono::steady_clock::time_point time_start, time;
        uint64_t baseline_memory_alloc = 0;
        uint64_t target_index_size = 0;
        bool log = false;
        uint16_t p = 0;
        pos_t roll_threshold = 0;

        char_t* T;
        pos_t n = 0;
        pos_t size_sss = 0;
        pos_t num_lpf = 0;
        pos_t len_lpf_phr = 0;
        pos_t num_fact = 0;
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

        struct lpf_pos_t {
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

        factorizer(char_t* input, pos_t input_size, parameters params)
            : log(params.log)
            , p(params.num_threads)
            , T(input)
            , n(input_size)
        { }

        template <
            quality_mode qual_mode,
            factorize_mode fact_mode = default_fact_mode,
            phrase_mode phr_mode = default_phr_mode,
            transform_mode transf_mode = default_transf_mode,
            template <typename> typename range_ds_t = default_range_ds_t,
            typename output_fnc_t>
        void factorize(output_fnc_t output)
        {
            static_assert(sizeof(char_t) == 1);

            if (p == 0) {
                p = omp_get_max_threads();
            }

            baseline_memory_alloc = malloc_count_current();
            malloc_count_reset_peak();
            omp_set_num_threads(p);

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

            if constexpr (qual_mode == exact) {
                static_assert(fact_mode != skip_phrases);
                std::string aprx_file_name = std::filesystem::temp_directory_path().string()
                    + "/aprx_" + random_alphanumeric_string(10);
                std::ofstream ofile_aprx(aprx_file_name);
                std::ostream_iterator<factor> ofile_aprx_it(ofile_aprx, "");
                compute_approximation<fact_mode, phr_mode>([&](factor f) { *ofile_aprx_it++ = f; });
                ofile_aprx.close();
                pos_t delta = std::min<pos_t>(n / num_fact, max_delta);
                pos_t max_num_samples = num_fact + n / delta;

                if (std::is_same_v<pos_t, uint32_t> ||
                    max_num_samples <= std::numeric_limits<uint32_t>::max()
                ) {
                    exact_factorizer<uint32_t, transf_mode, range_ds_t>(
                        T, n, LCE, aprx_file_name, delta, num_fact, p, log)
                        .transform_to_exact(output);
                } else {
                    exact_factorizer<uint64_t, transf_mode, range_ds_t>(
                        T, n, LCE, aprx_file_name, delta, num_fact, p, log)
                        .transform_to_exact(output);
                }

                std::filesystem::remove(aprx_file_name);
            } else {
                compute_approximation<fact_mode, phr_mode>(output);
            }

            if (log && fact_mode != skip_phrases) {
                uint64_t time_total = time_diff_ns(time_start, now());
                uint64_t mem_peak = malloc_count_peak() - baseline_memory_alloc;
                double comp_ratio = n / (double) num_fact;

                std::cout << "num. of factors: " << num_fact << std::endl;
                std::cout << "input length / num. of factors: " << comp_ratio << std::endl;
                std::cout << "total time: " << format_time(time_total) << std::endl;
                std::cout << "throughput: " << format_throughput(n, time_total) << std::endl;
                std::cout << "peak memory consumption: " << format_size(mem_peak) << std::endl;

                #ifdef LZ77_SSS_BENCH
                if (result_file_path != "") {
                    result_file
                        << " num_factors=" << num_fact
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
            phrase_mode phr_mode = default_phr_mode,
            typename output_fnc_t>
        void compute_approximation(output_fnc_t output)
        {
            LPF.resize(p);

            if constexpr (phr_mode == lpf_naive) {
                build_lce();
                build_LPF_naive();
            } else if constexpr (phr_mode == lpf_opt) {
                build_lce();
                build_LPF_opt();
            } else {
                if (log) std::cout << "reversing input" << std::flush;
                std::reverse(T, T + n);
                if (log) time = log_runtime(time);
                build_lce();
                build_LNF_all<phr_mode>();
                LCE = lce_t();
                if (log) std::cout << "reversing input" << std::flush;
                std::reverse(T, T + n);
                if (log) time = log_runtime(time);
                build_lce();
                build_LPF_all<phr_mode>();
            }

            LCE.delete_sa_s();

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
                double patt_len_guess = get_patt_len_guess(avg_gap_len, avg_lpf_phr_len, rel_len_gaps);

                if (log) {
                    log_phase("phrase_info", time_diff_ns(time, now()));
                    time = log_runtime(time);
                    std::cout << "|S| / (2n / tau) = " << size_sss / ((2.0 * n) / tau) << std::endl;
                    std::cout << "the density condition has " << (LCE.has_runs() ? "" : "not ") << "been applied" << std::endl;
                    std::cout << "num. of LPF phrases / SSS size = " << lpf_phr_per_sync << std::endl;
                    std::cout << "gaps length / input length: " << rel_len_gaps << std::endl;
                    std::cout << "num. of gaps / num. of LPF phrases: " << gaps_per_lpf_phr << std::endl;
                    std::cout << "avg. gap length: " << avg_gap_len << std::endl;
                    std::cout << "avg. LPF phrase length: " << avg_lpf_phr_len << std::endl;
                    std::cout << "pattern length guess: " << patt_len_guess << std::endl;
                    std::cout << "peak memory consumption: " << format_size(malloc_count_peak() - baseline_memory_alloc) << std::endl;
                    std::cout << "current memory consumption: " << format_size(malloc_count_current() - baseline_memory_alloc) << std::endl;
                    std::cout << "target index size: " << format_size(target_index_size) << std::endl;
                }

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

                if (fact_mode == greedy && !LCE.has_runs() && size_sss < 1.3 * ((2.0 * n) / tau) &&
                    n > min_par_input_size && rel_len_gaps > min_par_rel_gap_len && p > 1
                ) {
                    par_gap_idx = par_gap_idx_t(T, n, patt_lens, target_index_size, p);
                    if (log) std::cout << " (size: " << format_size(par_gap_idx.size_in_bytes()) << ")";
                } else {
                    gap_idx = gap_idx_t(T, n, patt_lens, target_index_size);
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
            return lce_l_64<pos_t>(T, i, j, max_lce);
        }

        static void greedy_phrase_selection(std::vector<lpf>& P);

        void build_lce();

        void build_PSV_NSV_S();

        void build_PGV_NGV_S();

        void build_LPF_naive();

        void build_LPF_opt();

        template <phrase_mode phr_mode>
        void build_LNF_all();
        
        template <phrase_mode phr_mode>
        void build_LPF_all();

        void get_phrase_info();

        inline factor longest_prev_occ(pos_t pos);

        template <bool first_block>
        inline factor longest_prev_occ_par(fp_arr_t& fps, pos_t pos, pos_t blk_end);

        template <bool first_block, typename next_lpf_t>
        void factorize_block(next_lpf_t next_lpf, pos_t blk_beg, pos_t blk_end);

        template <factorize_mode fact_mode, typename output_fnc_t>
        void factorize(output_fnc_t output);

        template <factorize_mode fact_mode, typename output_fnc_t, typename lpf_beg_t, typename next_lpf_t>
        void factorize_sequential(output_fnc_t output, lpf_beg_t lpf_beg, next_lpf_t next_lpf)
        {
            if constexpr (fact_mode == skip_phrases) {
                factorize_skip_gaps(output, lpf_beg, next_lpf);
            } else if constexpr (fact_mode == greedy) {
                factorize_greedy(output, lpf_beg, next_lpf);
            } else if constexpr (fact_mode == greedy_naive) {
                factorize_greedy_naive(output, lpf_beg, next_lpf);
            }
        }

        template <typename output_fnc_t, typename lpf_beg_t, typename next_lpf_t>
        void factorize_skip_gaps(output_fnc_t output, lpf_beg_t lpf_beg, next_lpf_t next_lpf);

        template <typename output_fnc_t, typename lpf_beg_t, typename next_lpf_t>
        void factorize_greedy_naive(output_fnc_t output, lpf_beg_t lpf_beg, next_lpf_t next_lpf);

        template <typename output_fnc_t, typename lpf_beg_t, typename next_lpf_t>
        void factorize_greedy(output_fnc_t output, lpf_beg_t lpf_beg, next_lpf_t next_lpf);

        template <typename output_fnc_t, typename lpf_beg_t, typename next_lpf_t>
        void factorize_greedy_parallel(output_fnc_t output, lpf_beg_t lpf_beg, next_lpf_t next_lpf);

        template <
            typename sidx_t,
            transform_mode transf_mode,
            template <typename> typename range_ds_t>
        class exact_factorizer {
        public:
            using sample_index_t = sample_index<pos_t, sidx_t, char_t, lce_t>;
            using point_t = typename range_ds_t<sidx_t>::point_t;
            using interval_t = sample_index_t::interval_t;
            using query_context_t = sample_index_t::query_ctx_t;

            std::chrono::steady_clock::time_point time_start, time;
            std::string aprx_file_name;
            std::string fact_file_name;
            bool log = false;
            uint16_t p = 0;

            char_t* T;
            const lce_t& LCE;

            pos_t n = 0;
            pos_t c = 0;
            pos_t delta = 0;
            pos_t& num_fact;

            struct sect_info_t {
                pos_t beg;
                sidx_t phr_idx;
            };

            std::vector<sect_info_t> par_sect;

            std::vector<pos_t> C;
            sample_index_t idx_C;
            std::vector<point_t> P;
            range_ds_t<sidx_t> R;
            std::vector<sidx_t> Pi;
            std::vector<sidx_t> Psi;

            inline pos_t LCE_R(pos_t i, pos_t j)
            {
                return LCE.lce(i, j);
            }

            exact_factorizer(char_t* T, pos_t n, const lce_t& LCE, std::string aprx_file_name,
                pos_t delta, pos_t& num_fact, uint16_t p, bool log
            ) : log(log), p(p), T(T), LCE(LCE), aprx_file_name(aprx_file_name), n(n), delta(delta), num_fact(num_fact) { }

            template <typename output_fnc_t>
            void transform_to_exact(output_fnc_t output)
            {
                if (log) {
                    time = now();
                    time_start = time;
                }

                build_c();
                build_idx_C();
                build_p();

                if constexpr (transf_mode != naive) {
                    build_pi_psi();
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
                    par_sect.resize(p + 1);
                    par_sect[p] = {.beg = n, .phr_idx = num_fact};
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

            void build_c();

            void build_idx_C();

            void build_pi_psi();

            void build_p();

            void insert_points(sidx_t& x_c, pos_t i);

            void find_close_sources(factor& f, pos_t i, pos_t e);

            inline void adjust_xc(sidx_t& gap_idx, pos_t pos);

            bool intersect(
                const interval_t& pa_c_iv, const interval_t& sa_c_iv,
                pos_t i, pos_t j, pos_t lce_l, pos_t lce_r, sidx_t& x_c, factor& f);

            template <typename output_fnc_t>
            void transform_to_exact_naive(output_fnc_t output);

            template <typename output_fnc_t>
            void transform_to_exact_without_samples(output_fnc_t output);

            void extend_right_with_samples(
                const interval_t& pa_c_iv,
                pos_t i, pos_t j, pos_t e, sidx_t& x_c, factor& f);

            template <typename output_fnc_t>
            void transform_to_exact_with_samples(output_fnc_t output);

            template <typename output_fnc_t>
            void combine_factorizations(output_fnc_t output);
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