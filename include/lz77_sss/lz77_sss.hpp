#pragma once

#include <vector>
#include <functional>
#include <filesystem>

#include <lce/lce_naive_wordwise_xor.hpp>
#include <lce/lce_sss.hpp>
#include <lce/lce_classic.hpp>

#include <lz77_sss/algorithms/lce_l.hpp>
#include <lz77_sss/misc/utils.hpp>
#include <lz77_sss/data_structures/rolling_hash_index_61.hpp>
#include <lz77_sss/data_structures/rolling_hash_index_107.hpp>
#include <lz77_sss/data_structures/sample_index/sample_index.hpp>
#include <lz77_sss/data_structures/dynamic_range/dynamic_square_grid.hpp>
#include <lz77_sss/data_structures/dynamic_range/semi_dynamic_square_grid.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_kd_tree.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_square_grid.hpp>
#include <lz77_sss/data_structures/static_weighted_range/static_weighted_striped_square.hpp>

enum phrase_mode {
    lpf_naive,
    lpf_all,
    lpf_all_external,
    lpf_lnf_all,
};

enum factorize_mode {
    greedy,
    greedy_naive,
    blockwise_all
};

enum transform_mode {
    naive,
    with_samples,
    without_samples,
};

template <typename pos_t = uint32_t>
class lz77_sss {
    public:
    static_assert(std::is_same_v<pos_t, uint32_t> || std::is_same_v<pos_t, uint64_t>);

    static constexpr uint64_t           default_tau              = 512;
    static constexpr phrase_mode        default_phr_mode         = lpf_all_external;
    static constexpr factorize_mode     default_fact_mode        = greedy;
    static constexpr transform_mode     default_transf_mode      = without_samples;
    template <typename sidx_t> using    default_range_ds_t       = decomposed_static_weighted_square_grid<sidx_t>;

    static constexpr uint8_t            num_patt_lens            = 5;
    static constexpr pos_t              delta_max                = 256;
    static constexpr pos_t              range_scan_threshold     = 4096;
    static constexpr uint64_t           max_buffered_lpf_phrases = 2048;

    struct factor {
        pos_t src;
        pos_t len;

        friend class lz77_sss;

        friend std::istream& operator>>(std::istream& in, const factor& f) {
            in.read((char*) &f, sizeof(factor));
            return in;
        }

        friend std::ostream& operator<<(std::ostream& out, const factor& f) {
            out.write((char*) &f, sizeof(factor));
            return out;
        }
    };

    template <
        factorize_mode                  fact_mode   = default_fact_mode,
        phrase_mode                     phr_mode    = default_phr_mode,
        uint64_t                        tau         = default_tau,
        std::output_iterator<factor>    out_it_t
    >
    static void factorize_approximate(std::string& input, out_it_t out_it, bool log = false) {
        factorizer<tau>(input, log).template factorize<approximate, fact_mode, phr_mode>(out_it);
    }

    template <
        factorize_mode  fact_mode   = default_fact_mode,
        phrase_mode     phr_mode    = default_phr_mode,
        uint64_t        tau         = default_tau
    >
    static void factorize_approximate(std::string& input, std::ofstream& out, bool log = false) {
        factorize_approximate<fact_mode, phr_mode, tau>(input, std::ostream_iterator<factor>(out, ""), log);
    }

    template <
        factorize_mode  fact_mode   = default_fact_mode,
        phrase_mode     phr_mode    = default_phr_mode,
        uint64_t        tau         = default_tau
    >
    static std::vector<factor> factorize_approximate(std::string& input, bool log = false) {
        std::vector<factor> factorization;
        factorize_approximate<fact_mode, phr_mode, tau>(input, std::back_insert_iterator(factorization), log);
        return factorization;
    }

    template <
        factorize_mode                  fact_mode       = default_fact_mode,
        phrase_mode                     phr_mode        = default_phr_mode,
        transform_mode                  transf_mode     = default_transf_mode,
        template <typename> typename    range_ds_t      = default_range_ds_t,
        uint64_t                        tau             = default_tau,
        std::output_iterator<factor>    out_it_t
    >
    static void factorize_exact(std::string& input, out_it_t out_it, bool log = false) {
        factorizer<tau>(input, log).template factorize<exact, fact_mode, phr_mode, transf_mode, range_ds_t>(out_it);
    }

    template <
        factorize_mode                  fact_mode       = default_fact_mode,
        phrase_mode                     phr_mode        = default_phr_mode,
        transform_mode                  transf_mode     = default_transf_mode,
        template <typename> typename    range_ds_t      = default_range_ds_t,
        uint64_t                        tau             = default_tau
    >
    static void factorize_exact(std::string& input, std::ofstream& out, bool log = false) {
        factorize_exact<fact_mode, phr_mode, transf_mode, range_ds_t, tau>(
            input, std::ostream_iterator<factor>(out, ""), log);
    }

    template <
        factorize_mode                  fact_mode       = default_fact_mode,
        phrase_mode                     phr_mode        = default_phr_mode,
        transform_mode                  transf_mode     = default_transf_mode,
        template <typename> typename    range_ds_t      = default_range_ds_t,
        uint64_t                        tau             = default_tau
    >
    static std::vector<factor> factorize_exact(std::string& input, bool log = false) {
        std::vector<factor> factorization;
        factorize_exact<fact_mode, phr_mode, transf_mode, range_ds_t, tau>(
            input, std::back_insert_iterator(factorization), log);
        return factorization;
    }
    
    template <std::input_iterator fact_it_t>
    static std::string decode(fact_it_t fact_it, pos_t output_size);

    static std::string decode(std::vector<factor>& factorization, pos_t output_size) {
        return decode(factorization.begin(), output_size);
    }

    static std::string decode(std::ifstream& in, pos_t output_size) {
        return decode(std::istream_iterator<factor>(in), output_size);
    }

    protected:

    struct lpf {
        pos_t beg;
        pos_t end;
        pos_t src;
        
        friend std::istream& operator>>(std::istream& in, const lpf& p) {
            in.read((char*) &p, sizeof(lpf));
            return in;
        }

        friend std::ostream& operator<<(std::ostream& out, const lpf& p) {
            out.write((char*) &p, sizeof(lpf));
            return out;
        }
    };

    enum lpf_mode {
        naive,
        all
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
        std::chrono::steady_clock::time_point time_start, time, time_end;
        uint64_t baseline_memory_alloc = 0;
        uint64_t target_index_size = 0;
        bool log = false;
        std::string lpf_file_name, sel_lpf_file_name;
        
        std::string& T;
        pos_t n = 0;
        uint32_t size_sss = 0;
        uint32_t num_lpf = 0;
        pos_t len_lpf_phr = 0;
        pos_t num_phr = 0;
        pos_t len_gaps = 0;
        pos_t num_gaps = 0;

        lce_t LCE;
        std::vector<lpf> LPF;
        std::array<pos_t, num_patt_lens> patt_lens;
        gap_idx_t gap_idx;

        std::vector<uint32_t> PSV_S;
        std::vector<uint32_t> NSV_S;
        std::vector<uint32_t> PGV_S;
        std::vector<uint32_t> NGV_S;

        factorizer(std::string& input, bool log) : log(log), T(input), n(input.size()) {}

        template <
            quality_mode                    qual_mode,
            factorize_mode                  fact_mode       = default_fact_mode,
            phrase_mode                     phr_mode        = default_phr_mode,
            transform_mode                  transf_mode     = default_transf_mode,
            template <typename> typename    range_ds_t      = default_range_ds_t,
            typename                        out_it_t
        >
        void factorize(out_it_t& out_it) {
            if (log) {
                time = now();
                time_start = time;
                baseline_memory_alloc = malloc_count_current();
                malloc_count_reset_peak();
            }

            omp_set_num_threads(1);

            if constexpr (qual_mode == exact) {
                std::string file_aprx_name = "aprx_" + random_alphanumeric_string(10);
                std::ofstream ofile_aprx(file_aprx_name);
                std::ostream_iterator<factor> ofile_aprx_it(ofile_aprx, "");
                compute_approximation<fact_mode, phr_mode>(ofile_aprx_it);
                ofile_aprx.close();
                std::ifstream ifile_aprx(file_aprx_name);
                std::istream_iterator<factor> ifile_aprx_it(ifile_aprx);
                pos_t delta = std::min<pos_t>(n / num_phr, delta_max);

                if constexpr (std::is_same_v<pos_t, uint32_t>) {
                    exact_factorizer<uint32_t, transf_mode, range_ds_t, out_it_t>(
                        T, LCE, delta, num_phr, log).transform_to_exact(ifile_aprx_it, out_it);
                } else {
                    pos_t max_num_samples = num_phr + n / delta;

                    if (max_num_samples <= std::numeric_limits<uint32_t>::max()) {
                        exact_factorizer<uint32_t, transf_mode, range_ds_t, out_it_t>(
                            T, LCE, delta, num_phr, log).transform_to_exact(ifile_aprx_it, out_it);
                    } else {
                        exact_factorizer<uint64_t, transf_mode, range_ds_t, out_it_t>(
                            T, LCE, delta, num_phr, log).transform_to_exact(ifile_aprx_it, out_it);
                    }
                }

                ifile_aprx.close();
                std::filesystem::remove(file_aprx_name);
            } else {
                compute_approximation<fact_mode, phr_mode>(out_it);
            }

            if (log) {
                uint64_t time_total = time_diff_ns(time_start, now());
                uint64_t mem_peak = malloc_count_peak() - baseline_memory_alloc;
                double comp_ratio = n / (double) num_phr;

                std::cout << "compression ratio: " << comp_ratio << std::endl;
                std::cout << "total time: " << format_time(time_total) << std::endl;
                std::cout << "throughput: " << format_throughput(n, time_total) << std::endl;
                std::cout << "peak memory consumption: " << format_size(mem_peak) << std::endl;

                #ifdef LZ77_SSS_BENCH
                if (result_file_path != "") {
                    std::ofstream result_file(
                        result_file_path, std::ofstream::app);
                    uint16_t transf_mode_int = qual_mode !=
                        exact ? 0 : (transf_mode + 1);

                    result_file << "RESULT"
                        << " text_name=" << text_name
                        << " n=" << n
                        << " alg=lz77_sss"
                        << " tau=" << tau
                        << " phr_mode=" << phr_mode
                        << " fact_mode=" << fact_mode
                        << " transf_mode=" << transf_mode_int
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
            factorize_mode  fact_mode   = default_fact_mode,
            phrase_mode     phr_mode    = default_phr_mode,
            typename        out_it_t
        >
        void compute_approximation(out_it_t& out_it) {
            if constexpr (phr_mode == lpf_naive) {
                build_lce();
                build_LPF_greedy();
            } else if constexpr (phr_mode == lpf_all) {
                build_lce();
                build_LPF_all([&](lpf&& p, pos_t){LPF.emplace_back(p);});
            } else if constexpr (phr_mode == lpf_all_external) {
                build_lce();
                build_LPF_all_external();
            } else {
                if (log) std::cout << "reversing T" << std::flush;
                std::reverse(T.begin(), T.end());
                if (log) time = log_runtime(time);
                build_lce();
                build_LNF_all([&](lpf&& p, pos_t){LPF.emplace_back(p);});
                LCE = lce_t();
                if (log) std::cout << "reversing T" << std::flush;
                std::reverse(T.begin(), T.end());
                if (log) time = log_runtime(time);
                build_lce();
                build_LPF_all([&](lpf&& p, pos_t){LPF.emplace_back(p);});
            }

            LCE.delete_ssa();

            if constexpr (phr_mode != lpf_naive) {
                if (log) {
                    std::cout << "greedily selecting LPF phrases" << std::flush;
                }

                if constexpr (phr_mode == lpf_all_external) {
                    greedy_phrase_selection_external();
                } else {
                    greedy_phrase_selection(LPF);
                }

                if (log) {
                    time = log_runtime(time);
                }
            }

            if (log) std::cout << "inspecting LPF" << std::flush;

            if constexpr (phr_mode == lpf_all_external) {
                get_phrase_info_external();
                std::ofstream lpf_ofile(sel_lpf_file_name, std::ios::app);
                lpf_ofile << lpf {.beg = n, .end = n + 1};
                lpf_ofile.close();
            } else {
                get_phrase_info();
                LPF.emplace_back(lpf {.beg = n, .end = n + 1});
                LPF.shrink_to_fit();
            }

            double lpf_phr_per_sync = num_lpf / (double) size_sss;
            double rel_len_gaps = len_gaps / (double) n;
            double gaps_per_lpf_phr = num_gaps / (double) num_lpf;
            double avg_gap_len = len_gaps / (double) num_gaps;
            double avg_lpf_phr_len = len_lpf_phr / (double) num_lpf;
            
            target_index_size = std::max<uint64_t>({
                malloc_count_peak(),
                baseline_memory_alloc + ((n / 3.0) * rel_len_gaps),
                baseline_memory_alloc + (n / 10.0)
            }) - malloc_count_current();

            if (log) {
                time = log_runtime(time);
                std::cout << "num. of LPF phrases / SSS size = " << lpf_phr_per_sync << std::endl;
                std::cout << "gaps length / input length: " << rel_len_gaps << std::endl;
                std::cout << "num. of gaps / num. of LPF phrases: " << gaps_per_lpf_phr << std::endl;
                std::cout << "avg. gap length: " << avg_gap_len << std::endl;
                std::cout << "avg. LPF phrase length: " << avg_lpf_phr_len << std::endl;
                std::cout << "peak memory consumption: "
                    << format_size(malloc_count_peak() - baseline_memory_alloc) << std::endl;
                std::cout << "current memory consumption: "
                    << format_size(malloc_count_current() - baseline_memory_alloc) << std::endl;
                std::cout << "target index size: " << format_size(target_index_size) << std::endl;
            }

            double patt_len_guess = std::min<double>({
                avg_gap_len, avg_lpf_phr_len,
                8.0 * std::pow(128, 1.0 - rel_len_gaps)});
            
                 if (patt_len_guess <= 6)   {patt_lens = {2,3, 4, 5, 6};}
            else if (patt_len_guess <= 8)   {patt_lens = {2,3, 4, 6, 8};}
            else if (patt_len_guess <= 12)  {patt_lens = {2,3, 4, 8,12};}
            else if (patt_len_guess <= 16)  {patt_lens = {2,4, 6, 9,16};}
            else if (patt_len_guess <= 32)  {patt_lens = {2,4, 6,10,20};}
            else if (patt_len_guess <= 64)  {patt_lens = {2,4, 7,12,28};}
            else if (patt_len_guess <= 128) {patt_lens = {2,4, 8,16,36};}
            else if (patt_len_guess <= 256) {patt_lens = {2,5,10,20,42};}
            else if (patt_len_guess <= 512) {patt_lens = {2,6,12,24,48};}
            else                            {patt_lens = {2,8,16,32,64};}
            
            if (log) {
                std::cout << "pattern lengths for the rolling hash index: ";
                for (pos_t i = 0; i < num_patt_lens - 1; i++) std::cout << patt_lens[i] << ", ";
                std::cout << patt_lens[num_patt_lens - 1] << std::endl;
                std::cout << "initializing rolling hash index" << std::flush;
            }

            gap_idx = gap_idx_t(T.data(), n, patt_lens, target_index_size);

            if (log) {
                std::cout << " (size: " << format_size(gap_idx.size_in_bytes()) << ")";
                time = log_runtime(time);
            }

            if constexpr (phr_mode == lpf_all_external) {
                std::ifstream lpf_ifile(sel_lpf_file_name);
                std::istream_iterator<lpf> lpf_it(lpf_ifile);
                factorize<fact_mode>(out_it, [&](){return *lpf_it++;});
                std::filesystem::remove(sel_lpf_file_name);
            } else {
                uint32_t p = 0;
                factorize<fact_mode>(out_it, [&](){return LPF[p++];});
            }

            if (log) {
                #ifndef NDEBUG
                std::cout << "rate of initialized values"
                    << "in the rolling hash index: "
                    << gap_idx.rate_init()
                    << std::endl;
                #endif
            }

            LPF.clear();
            LPF.shrink_to_fit();
            gap_idx = gap_idx_t();
        }

        inline pos_t LCE_R(pos_t i, pos_t j) {
            return LCE.lce(i, j);
        }

        inline pos_t LCE_L(pos_t i, pos_t j, pos_t max_lce) {
            return lce_l_128<pos_t>(T.data(), i, j, max_lce);
        }

        static void greedy_phrase_selection(std::vector<lpf>& P);

        void greedy_phrase_selection_external();

        void build_lce();

        void build_PSV_NSV_S();
        
        void build_PGV_NGV_S();

        void build_LPF_greedy();

        void build_LNF_greedy();

        void build_LPF_all(std::function<void(lpf&&,pos_t)> lpf_it);

        void build_LNF_all(std::function<void(lpf&&,pos_t)> lpf_it);

        void build_LPF_all_external();

        void get_phrase_info();

        void get_phrase_info_external();

        inline factor longest_prev_occ(pos_t pos);

        template <factorize_mode fact_mode, typename out_it_t>
        void factorize(out_it_t& out_it, std::function<lpf()> lpf_it) {
            if constexpr (fact_mode == greedy)        {factorize_greedy<       >(out_it, lpf_it);} else
            if constexpr (fact_mode == greedy_naive)  {factorize_greedy_naive< >(out_it, lpf_it);} else
            if constexpr (fact_mode == blockwise_all) {factorize_blockwise_all<>(out_it, lpf_it);}
        }

        template <typename out_it_t>
        void factorize_greedy_naive(out_it_t& out_it, std::function<lpf()> lpf_it);

        template <typename out_it_t>
        void factorize_greedy(out_it_t& out_it, std::function<lpf()> lpf_it);

        template <typename out_it_t>
        void factorize_blockwise_all(out_it_t& out_it, std::function<lpf()> lpf_it);

        template <
            typename sidx_t,
            transform_mode transf_mode,
            template <typename> typename range_ds_t,
            typename out_it_t
        >
        class exact_factorizer {
            public:

            using sample_index_t = sample_index<pos_t, sidx_t, lce_t>;
            using point_t = typename range_ds_t<sidx_t>::point_t;
            using sxa_interval_t = sample_index<pos_t, sidx_t, lce_t>::sxa_interval_t;
            using query_context_t = sample_index<pos_t, sidx_t, lce_t>::query_context;
            using time_point_t = std::chrono::steady_clock::time_point;

            std::chrono::steady_clock::time_point time_start, time, time_end;
            bool log = false;

            const std::string& T;
            const lce_t& LCE;

            pos_t n = 0;
            pos_t c = 0;
            pos_t delta = 0;
            pos_t& num_phr;

            std::vector<pos_t> C;
            sample_index_t idx_C;
            std::vector<point_t> P;
            range_ds_t<sidx_t> R;
            std::vector<sidx_t> PS;
            std::vector<sidx_t> SP;

            uint64_t time_extend_left = 0;
            uint64_t time_extend_right = 0;
            uint64_t time_range_queries = 0;
            uint64_t time_close_sources = 0;
            uint64_t time_insert_points = 0;
            uint64_t num_range_queries = 0;
            uint64_t spa_range_sum = 0;
            uint64_t ssa_range_sum = 0;

            inline pos_t LCE_R(pos_t i, pos_t j) {
                return LCE.lce(i, j);
            }

            exact_factorizer(std::string& T, const lce_t& LCE, pos_t delta, pos_t& num_phr, bool log)
                : log(log), T(T), LCE(LCE), n(T.size()), delta(delta), num_phr(num_phr) {}

            void transform_to_exact(std::istream_iterator<factor>& ifile_aprx_it, out_it_t& out_it) {
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
                    R = range_ds_t<sidx_t>(T, C, P);
                } else {
                    R = range_ds_t<sidx_t>(P, c);
                }

                if constexpr (range_ds_t<sidx_t>::is_static()) {
                    P.clear();
                    P.shrink_to_fit();
                }

                if (log) {
                    std::cout << " (" << format_size(R.size_in_bytes()) << ")";
                    time = log_runtime(time);
                }

                if constexpr (transf_mode == naive) {
                    transform_to_exact_naive(out_it);
                } else if constexpr (transf_mode == with_samples) {
                    transform_to_exact_with_samples(out_it);
                } else if constexpr (transf_mode == without_samples) {
                    transform_to_exact_without_samples(out_it);
                }

                if (log) {
                    std::cout << "time for extend left: " << format_time(time_extend_left) << std::endl;
                    std::cout << "time for extend right: " << format_time(time_extend_right) << std::endl;
                    std::cout << "time for range queries: " << format_time(time_range_queries) << std::endl;
                    std::cout << "avg. SPA query range: " << spa_range_sum / (double) num_range_queries << std::endl;
                    std::cout << "avg. SSA query range: " << ssa_range_sum / (double) num_range_queries << std::endl;

                    if constexpr (range_ds_t<sidx_t>::is_dynamic()) {
                        std::cout << "time for finding close sources: " << format_time(time_close_sources) << std::endl;
                        std::cout << "time for inserting points: " << format_time(time_insert_points) << std::endl;
                        std::cout << "final size of " << range_ds_t<sidx_t>::name()
                            << ": " << format_size(R.size_in_bytes()) << std::endl;
                    }
                }
            }

            void build_c(std::istream_iterator<factor>& ifile_aprx_it);

            void build_idx_C();

            void build_ps_sp();

            void build_p();

            void insert_points(sidx_t& x_c, pos_t i);

            void handle_close_sources(factor& f, pos_t i);

            inline void adjust_xc(sidx_t& gap_idx, pos_t pos);

            bool intersect(
                const sxa_interval_t& spa_iv, const sxa_interval_t& ssa_iv,
                pos_t i, pos_t j, pos_t lce_l,pos_t lce_r, sidx_t& x_c, factor& f
            );

            void transform_to_exact_naive(out_it_t& out_it);

            void transform_to_exact_without_samples(out_it_t& out_it);

            void extend_right_with_samples(
                const sxa_interval_t& spa_iv,
                pos_t i, pos_t j, sidx_t& x_c, factor& f
            );

            void transform_to_exact_with_samples(out_it_t& out_it);
        };
    };
};

#include "algorithms/common.cpp"

#include "algorithms/approximate/common.cpp"
#include "algorithms/approximate/lpf_lnf/nxv_pxv.cpp"
#include "algorithms/approximate/lpf_lnf/greedy.cpp"
#include "algorithms/approximate/lpf_lnf/all.cpp"
#include "algorithms/approximate/lpf_lnf/external.cpp"
#include "algorithms/approximate/factorize/common.cpp"
#include "algorithms/approximate/factorize/greedy.cpp"
#include "algorithms/approximate/factorize/greedy_naive.cpp"
#include "algorithms/approximate/factorize/blockwise_all.cpp"

#include "algorithms/transform_to_exact/common.cpp"
#include "algorithms/transform_to_exact/naive.cpp"
#include "algorithms/transform_to_exact/with_samples.cpp"
#include "algorithms/transform_to_exact/without_samples.cpp"