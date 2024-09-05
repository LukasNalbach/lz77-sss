#include <fstream>

std::string text_name;
std::string result_file_path;
#define LZ77_SSS_BENCH 1

#include <lz77_sss/lz77_sss.hpp>
#include <lz77/lpf_factorizer.hpp>
#include <lz77/gzip9_factorizer.hpp>

std::string T;
uint64_t n;
std::ifstream input_file;

template <typename pos_t>
void check_correctness(std::string file_name) {
    std::ifstream fact_sss_file(file_name);
    std::string T_rev = lz77_sss<pos_t>::decode(fact_sss_file, n);
    fact_sss_file.close();
    bool correct = true;

    if (T_rev.size() != n) {
        correct = false;
    } else {
        for (pos_t i = 0; i < n; i++) {
            if (T_rev[i] != T[i]) {
                correct = false;
                break;
            }
        }
    }

    std::cout << "the factorization is " <<
        (correct ? "" : "not ") << "correct" << std::endl;
}

template <
    factorize_mode fact_mode,
    phrase_mode phr_mode
> void run_sss_approximate(std::string file_name) {
    std::ofstream fact_sss_file(file_name);

    if (n <= std::numeric_limits<uint32_t>::max()) {
        lz77_sss<uint32_t>::factorize_approximate<
            fact_mode, phr_mode>(T, fact_sss_file, true);
        fact_sss_file.close();
        check_correctness<uint32_t>(file_name);
    } else {
        lz77_sss<uint64_t>::factorize_approximate<
            fact_mode, phr_mode>(T, fact_sss_file, true);
        fact_sss_file.close();
        check_correctness<uint64_t>(file_name);
    }
}

template <
    factorize_mode fact_mode,
    phrase_mode phr_mode,
    transform_mode transf_mode,
    template <typename> typename range_ds_t
>
void run_sss_exact(std::string file_name) {
    std::ofstream fact_sss_file(file_name);
    
    if (n <= std::numeric_limits<uint32_t>::max()) {
        lz77_sss<uint32_t>::factorize_exact<
            fact_mode, phr_mode, transf_mode, range_ds_t>(
            T, fact_sss_file, true);
        fact_sss_file.close();
        check_correctness<uint32_t>(file_name);
    } else {
        lz77_sss<uint64_t>::factorize_exact<
            fact_mode, phr_mode, transf_mode, range_ds_t>(
            T, fact_sss_file, true);
        fact_sss_file.close();
        check_correctness<uint64_t>(file_name);
    }
}

void log_algorithm(
    std::string alg_name, uint64_t time,
    uint64_t mem_peak, uint64_t num_factors
) {
    double comp_ratio = n / (double) num_factors;

    std::cout << ", in ~ " << format_time(time) << std::endl;
    std::cout << "throughput: "
        << format_throughput(n, time) << std::endl;
    std::cout << "peak memory consumption: "
        << format_size(mem_peak) << std::endl;
    std::cout << "compression ratio: " << comp_ratio << std::endl;

    if (result_file_path != "") {
        std::ofstream result_file(
            result_file_path, std::ofstream::app);

        result_file << "RESULT"
            << " text_name=" << text_name
            << " n=" << n
            << " alg=" << alg_name
            << " num_factors=" << num_factors
            << " comp_ratio=" << comp_ratio
            << " time=" << time
            << " throughput=" << throughput(n, time)
            << " mem_peak=" << mem_peak << std::endl;
    }
}

int main(int argc, char** argv) {
    if (!(argc == 2 || argc == 3)) {
        std::cout << "usage: lz77_sss_bench <file> <result_file>";
        exit(-1);
    }

    std::string file_path = argv[1];
    input_file.open(file_path);
    
    if (!input_file.good()) {
        std::cout << "error: could not read <file>";
        exit(-1);
    }

    if (argc == 3) {
        text_name = file_path.substr(
            file_path.find_last_of("/\\") + 1);
        result_file_path = argv[2];
    }

    input_file.seekg(0, std::ios::end);
    n = input_file.tellg();
    input_file.seekg(0, std::ios::beg);
    auto t0 = now();
    std::cout << "reading T (" << format_size(n) << ")" << std::flush;
    no_init_resize(T, n);
    read_from_file(input_file, T.data(), n);
    input_file.close();
    log_runtime(t0);

    std::cout << std::endl << "running naive LZ77 SSS 3-approximation:" << std::endl;
    run_sss_approximate<greedy_naive, lpf_naive>("fact_sss_aprx");
    std::filesystem::remove("fact_sss_aprx");

    std::cout << std::endl << "running LZ77 SSS 3-approximation:" << std::endl;
    run_sss_approximate<greedy_optimized, lpf_optimal>("fact_sss_aprx");
    std::filesystem::remove("fact_sss_aprx");

    std::cout << std::endl << "running LZ77 SSS 1.5-approximation:" << std::endl;
    run_sss_approximate<greedy_optimized, lpf_lnf_optimal>("fact_sss_aprx");
    std::filesystem::remove("fact_sss_aprx");
    
    std::cout << std::endl << "running naive LZ77 SSS exact algorithm:" << std::endl;
    run_sss_exact<greedy_optimized, lpf_optimal, naive,
        semi_dynamic_square_grid>("fact_sss_exact");
    std::filesystem::remove("fact_sss_exact");

    std::cout << std::endl << "running LZ77 SSS exact algorithm (with samples):" << std::endl;
    run_sss_exact<greedy_optimized, lpf_optimal, optimized_with_samples,
        semi_dynamic_square_grid>("fact_sss_exact");
    std::filesystem::remove("fact_sss_exact");
    
    std::cout << std::endl << "running LZ77 SSS exact algorithm (without samples):" << std::endl;
    run_sss_exact<greedy_optimized, lpf_optimal, optimized_without_samples,
        semi_dynamic_square_grid>("fact_sss_exact");
    std::filesystem::remove("fact_sss_exact");

    std::cout << std::endl << "running LZ77 LPF algorithm" << std::flush;
    uint64_t baseline_memory_alloc = malloc_count_current();
    std::ofstream file_lpf("fact_lpf");
    malloc_count_reset_peak();
    auto t1 = now();
    lz77::LPFFactorizer().factorize(T.begin(), T.end(),
        std::ostream_iterator<lz77::Factor>(file_lpf, ""));
    auto t2 = now();
    uint64_t num_factors = file_lpf.tellp() / sizeof(lz77::Factor);
    uint64_t time = time_diff_ns(t1, t2);
    uint64_t mem_peak = malloc_count_peak() - baseline_memory_alloc;
    file_lpf.close();
    std::filesystem::remove("fact_lpf");
    log_algorithm("lpf", time, mem_peak, num_factors);

    std::cout << std::endl << "running gzip -9" << std::flush;
    std::ofstream file_gz9("fact_gz9");
    malloc_count_reset_peak();
    auto t3 = now();
    lz77::Gzip9Factorizer().factorize(T.begin(), T.end(),
        std::ostream_iterator<lz77::Factor>(file_gz9, ""));
    auto t4 = now();
    num_factors = file_gz9.tellp() / sizeof(lz77::Factor);
    time = time_diff_ns(t3, t4);
    mem_peak = malloc_count_peak() - baseline_memory_alloc;
    file_gz9.close();
    std::filesystem::remove("fact_gz9");
    log_algorithm("gz9", time, mem_peak, num_factors);
}