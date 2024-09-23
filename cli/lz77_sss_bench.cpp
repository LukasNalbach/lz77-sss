#include <fstream>

std::string text_name;
std::string result_file_path;
#define LZ77_SSS_BENCH 1

#include <lz77_sss/lz77_sss.hpp>
#include <lz77/lpf_factorizer.hpp>
#include <lz77/gzip9_factorizer.hpp>

uint16_t max_threads = 0;
std::string T;
uint64_t n;

template <
    factorize_mode fact_mode,
    phrase_mode phr_mode
> void run_sss_approximate(std::string file_name, uint16_t max_threads) {

    for (uint16_t num_threads = 1; num_threads < max_threads; num_threads *= 2) {
        std::ofstream fact_sss_file(file_name);

        if (n <= std::numeric_limits<uint32_t>::max()) {
            lz77_sss<uint32_t>::factorize_approximate<
                fact_mode, phr_mode>(T, fact_sss_file,
                {.num_threads = num_threads, .log = true});
        } else {
            lz77_sss<uint64_t>::factorize_approximate<
                fact_mode, phr_mode>(T, fact_sss_file,
                {.num_threads = num_threads, .log = true});
        }
    }
}

template <
    factorize_mode fact_mode,
    phrase_mode phr_mode,
    transform_mode transf_mode,
    template <typename> typename range_ds_t
>
void run_sss_exact(std::string file_name, uint16_t max_threads) {

    for (uint16_t num_threads = 1; num_threads < max_threads; num_threads *= 2) {
        std::ofstream fact_sss_file(file_name);

        if (n <= std::numeric_limits<uint32_t>::max()) {
            lz77_sss<uint32_t>::factorize_exact<
                fact_mode, phr_mode, transf_mode, range_ds_t>(
                T, fact_sss_file, {.num_threads = num_threads, .log = true});
        } else {
            lz77_sss<uint64_t>::factorize_exact<
                fact_mode, phr_mode, transf_mode, range_ds_t>(
                T, fact_sss_file, {.num_threads = num_threads, .log = true});
        }
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
    if (!(2 <= argc && argc <= 4)) {
        std::cout << "usage: lz77_sss_bench <file> <max_threads> <result_file>";
        std::cout << "       the last two parameters are optional";
        exit(-1);
    }

    std::string file_path = argv[1];
    std::ifstream input_file(file_path);
    
    if (!input_file.good()) {
        std::cout << "error: could not read <file>";
        exit(-1);
    }

    if (argc >= 3) {
        max_threads = atoi(argv[2]);

        if (max_threads == 0 || max_threads > omp_get_max_threads()) {
            std::cout << "error: invalid number of threads";
            exit(-1);
        }
    }

    if (argc == 4) {
        text_name = file_path.substr(
            file_path.find_last_of("/\\") + 1);
        result_file_path = argv[3];
    }

    input_file.seekg(0, std::ios::end);
    n = input_file.tellg();
    input_file.seekg(0, std::ios::beg);
    auto t0 = now();
    std::cout << "reading T (" << format_size(n) << ")" << std::flush;
    no_init_resize_with_exess(T, n, 4 * lz77_sss<>::default_tau);
    read_from_file(input_file, T.data(), n);
    input_file.close();
    log_runtime(t0);

    std::cout << std::endl << "running naive LZ77 SSS 3-approximation:" << std::endl;
    run_sss_approximate<greedy_naive, lpf_naive>("fact_sss_aprx", 1);
    std::filesystem::remove("fact_sss_aprx");

    std::cout << std::endl << "running LZ77 SSS 3-approximation:" << std::endl;
    run_sss_approximate<greedy, lpf_all_external>("fact_sss_aprx", max_threads);
    std::filesystem::remove("fact_sss_aprx");
    
    std::cout << std::endl << "running LZ77 SSS 1.5-approximation:" << std::endl;
    run_sss_approximate<greedy, lpf_lnf_all>("fact_sss_aprx", max_threads);
    std::filesystem::remove("fact_sss_aprx");
    
    std::cout << std::endl << "running naive LZ77 SSS exact algorithm:" << std::endl;
    run_sss_exact<greedy, lpf_all_external, naive,
        static_weighted_square_grid>("fact_sss_exact", 1);
    std::filesystem::remove("fact_sss_exact");
    
    std::cout << std::endl << "running LZ77 SSS exact algorithm (with samples):" << std::endl;
    run_sss_exact<greedy, lpf_all_external, with_samples,
        decomposed_static_weighted_square_grid>("fact_sss_exact", max_threads);
    std::filesystem::remove("fact_sss_exact");

    std::cout << std::endl << "running LZ77 SSS exact algorithm (without samples):" << std::endl;
    std::cout << std::endl << "range data structure: decomposed static weighted square grid" << std::endl;
    run_sss_exact<greedy, lpf_all_external, without_samples,
        decomposed_static_weighted_square_grid>("fact_sss_exact", max_threads);
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

    baseline_memory_alloc = malloc_count_current();
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