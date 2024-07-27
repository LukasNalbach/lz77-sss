#pragma once

#include <vector>
#include <functional>

//#include <absl/container/btree_set.h>
#include <lce/lce_sss.hpp>
#include <lce/lce_classic.hpp>

#include <lz77/lpf_factorizer.hpp>

#include <lz77_sss_approx/misc/utils.hpp>
#include <lz77_sss_approx/data_structures/rolling_hash_index_107.hpp>
#include <lz77_sss_approx/data_structures/rolling_hash_index_107.hpp>

enum phrase_mode {
    lpf_naive,
    lpf_optimal,
    lpf_lnf_optimal
};

enum factozize_mode {
    greedy,
    blockwise_optimal,
    hybrid
};

template <typename pos_t = uint32_t>
class lz77_sss_approx {
    public:
    static_assert(std::is_same_v<pos_t, uint32_t> || std::is_same_v<pos_t, uint64_t>);

    static constexpr uint64_t default_tau = 512;
    static constexpr phrase_mode default_phr_mode = lpf_optimal;
    static constexpr factozize_mode default_fact_mode = greedy;

    struct factor {
        pos_t src;
        pos_t len;

        friend class lz77_sss_approx;

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
        factozize_mode fact_mode = default_fact_mode,
        phrase_mode phr_mode = default_phr_mode,
        uint64_t tau = default_tau,
        std::output_iterator<factor> out_it_t
    > static void factorize(std::string& input, out_it_t out_it, bool log = false) {
        factorizer<fact_mode, phr_mode, tau, out_it_t>(input, out_it, log);
    }

    template <
        factozize_mode fact_mode = default_fact_mode,
        phrase_mode phr_mode = default_phr_mode,
        uint64_t tau = default_tau
    > static void factorize(std::string& input, std::ofstream& out, bool log = false) {
        factorize<fact_mode, phr_mode>(input, std::ostream_iterator<factor>(out, ""), log);
    }

    template <
        factozize_mode fact_mode = default_fact_mode,
        phrase_mode phr_mode = default_phr_mode,
        uint64_t tau = default_tau
    > static std::vector<factor> factorize(std::string& input, bool log = false) {
        std::vector<factor> factorization;
        factorize<fact_mode, phr_mode>(input, std::back_insert_iterator(factorization), log);
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

    enum lpf_mode {
        naive,
        all_extend_left
    };

    struct lpf {
        pos_t beg;
        pos_t end;
        pos_t src;
    };

    lz77_sss_approx() = delete;
    lz77_sss_approx(lz77_sss_approx&& other) = delete;
    lz77_sss_approx(const lz77_sss_approx& other) = delete;
    lz77_sss_approx& operator=(lz77_sss_approx&& other) = delete;
    lz77_sss_approx& operator=(const lz77_sss_approx& other) = delete;

    template <
        factozize_mode fact_mode,
        phrase_mode phr_mode,
        uint64_t tau,
        typename out_it_t
    > class factorizer {
        public:

        factorizer() = delete;
        factorizer(factorizer&& other) = delete;
        factorizer(const factorizer& other) = delete;
        factorizer& operator=(factorizer&& other) = delete;
        factorizer& operator=(const factorizer& other) = delete;

        using lce_t = alx::lce::lce_sss<char, tau, pos_t, false>;

        out_it_t& out_it;
        bool log = false;
        std::chrono::steady_clock::time_point time_start,time,time_end;
        uint64_t baseline_memory_alloc = 0;
        uint64_t target_index_size = 0;
        
        std::string& T;
        pos_t n = 0;
        uint8_t sigma = 0;
        uint8_t word_width = 0;
        uint32_t size_sss = 0;
        uint32_t num_lpf = 0;
        pos_t len_lpf_phr = 0;
        pos_t num_phr = 0;
        pos_t len_gaps = 0;
        pos_t num_gaps = 0;

        lce_t LCE;
        std::vector<pos_t> char_freq;
        std::vector<uint8_t> map_ext;
        std::vector<lpf> LPF;

        factorizer(std::string& input, out_it_t& out_it, bool log)
         : out_it(out_it), log(log), T(input), n(input.size())
        {
            if (log) {
                time = now();
                time_start = time;
                baseline_memory_alloc = malloc_count_current();
            }

            omp_set_num_threads(1);

            //compute_char_freq();
            //map_to_effective_alphabet();

            if constexpr (phr_mode == lpf_naive) {
                build_lce();
                build_LPF_S<naive>();
            } else if constexpr (phr_mode == lpf_optimal) {
                build_lce();
                build_LPF_S<all_extend_left>();
            } else if constexpr (phr_mode == lpf_lnf_optimal) {
                std::reverse(T.begin(),T.end());
                build_lce();
                build_LNF_S<all_extend_left>();
                LCE = lce_t();
                std::reverse(T.begin(),T.end());
                build_lce();
                build_LPF_S<all_extend_left>();
            }

            LCE.delete_ssa();

            if constexpr (phr_mode != lpf_naive) {
                if (log) {
                    std::cout << "greedily selecting LPF phrases" << std::flush;
                }

                greedy_phrase_selection(LPF);
                LPF.shrink_to_fit();

                if (log) {
                    time = log_runtime(time);
                }
            }

            get_phrase_info();
            
            double avg_gap_len = len_gaps / (double) num_gaps;
            double avg_lpf_phr_len = len_lpf_phr / (double) num_lpf;

            if (log) {
                std::cout << "num. of LPF phrases / SSS size = " << num_lpf / (double) size_sss << std::endl;
                std::cout << "gaps length / input length: " << len_gaps / (double) n << std::endl;
                std::cout << "num. of gaps / num. of LPF phrases: " << num_gaps / (double) num_lpf << std::endl;
                std::cout << "avg. gap length: " << avg_gap_len << std::endl;
                std::cout << "avg. LPF phrase length: " << avg_lpf_phr_len << std::endl;
            }
            
            target_index_size = std::max<uint64_t>(malloc_count_peak(), n / 3 + baseline_memory_alloc) - malloc_count_current();
            double patt_len_guess = std::min<double>(avg_gap_len, avg_lpf_phr_len);

                 if (patt_len_guess <= 12)  {factorize<fact_mode,2,3, 4, 8,12>();}
            else if (patt_len_guess <= 16)  {factorize<fact_mode,2,4, 6, 9,16>();}
            else if (patt_len_guess <= 32)  {factorize<fact_mode,2,4, 6,10,20>();}
            else if (patt_len_guess <= 64)  {factorize<fact_mode,2,4, 7,12,28>();}
            else if (patt_len_guess <= 128) {factorize<fact_mode,2,4, 8,16,36>();}
            else if (patt_len_guess <= 256) {factorize<fact_mode,2,5,10,20,42>();}
            else if (patt_len_guess <= 512) {factorize<fact_mode,2,6,12,24,48>();}
            else                            {factorize<fact_mode,2,8,16,32,64>();}

            //remap_to_original_alphabet();

            if (log) {
                std::cout << "compression ratio: " << n / (double) num_phr << std::endl;
                std::cout << "total time: " << format_time(time_diff_ns(time_start, time)) << std::endl;
                std::cout << "peak memory consumption: " << format_size(malloc_count_peak() - baseline_memory_alloc) << std::endl;
            }
        }

        inline pos_t LCE_R(const pos_t i, const pos_t j) {
            return LCE.lce(i,j);
        }

        pos_t LCE_L(pos_t i, pos_t j, pos_t max_lce = 32768);

        static void greedy_phrase_selection(std::vector<lpf>& P);

        void compute_char_freq();

        void map_to_effective_alphabet();

        void remap_to_original_alphabet();

        void build_lce();

        template <lpf_mode mode> void build_LPF_S();

        template <lpf_mode mode> void build_LNF_S();

        void get_phrase_info();

        template <pos_t... patt_lens>
        inline factor longest_prev_occ(rolling_hash_index_107<pos_t, patt_lens...>& idx, pos_t pos, pos_t len_max);

        template <pos_t... patt_lens>
        void factorize() {
            if constexpr (fact_mode == greedy)            {factorize_greedy            <patt_lens...>();} else
            if constexpr (fact_mode == blockwise_optimal) {factorize_blockwise_optimal <patt_lens...>();} else
            if constexpr (fact_mode == hybrid)            {factorize_hybrid            <patt_lens...>();} // WIP
        }

        template <pos_t... patt_lens>
        void factorize_greedy();

        template <pos_t... patt_lens>
        void factorize_blockwise_optimal();

        template <pos_t... patt_lens>
        inline lpf lpf_extend_left(rolling_hash_index_107<pos_t, patt_lens...>& idx, pos_t pos, pos_t max_lce_l);

        template <pos_t... patt_lens>
        void factorize_hybrid();
    };
};

#include "algorithms/common.cpp"
#include "algorithms/lpf_lnf.cpp"

#include "algorithms/factorize/blockwise_optimal.cpp"
#include "algorithms/factorize/greedy.cpp"
#include "algorithms/factorize/hybrid.cpp"