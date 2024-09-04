#include <fstream>

#include <lz77_sss/lz77_sss.hpp>
#include <lz77/lpf_factorizer.hpp>
#include <lz77/gzip9_factorizer.hpp>

std::string T;
uint64_t n;

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

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "usage: lz77_sss_bench <file>";
        exit(-1);
    }

    std::ifstream input_file;
    input_file.open(argv[1]);
    
    if (!input_file.good()) {
        std::cout << "error: could not read <file>";
        exit(-1);
    }

    input_file.seekg(0, std::ios::end);
    n = input_file.tellg();
    input_file.seekg(0, std::ios::beg);
    auto time = now();
    std::cout << "reading T (" << format_size(n) << ")" << std::flush;
    no_init_resize(T, n);
    read_from_file(input_file, T.data(), n);
    input_file.close();
    time = log_runtime(time);

    std::cout << std::endl << "running LZ77 SSS approximation:" << std::endl;
    run_sss_approximate<greedy_skip_phrases, lpf_optimal>("fact_sss_aprx");

    std::cout << std::endl << "running LZ77 SSS exact algorithm:" << std::endl;
    run_sss_exact<greedy_skip_phrases, lpf_optimal, optimized,
        semi_dynamic_square_grid>("fact_sss_exact");

    std::cout << std::endl << "running LZ77 LPF algorithm" << std::flush;
    uint64_t baseline_memory_alloc = malloc_count_current();
    std::ofstream file_exact("fact_exact");
    malloc_count_reset_peak();
    auto t1 = now();
    lz77::LPFFactorizer().factorize(T.begin(), T.end(),
        std::ostream_iterator<lz77::Factor>(file_exact, ""));
    auto t2 = now();
    uint64_t z = file_exact.tellp() / sizeof(lz77::Factor);
    file_exact.close();
    log_runtime(t1, t2);
    std::cout << "throughput: "
        << format_throughput(n, time_diff_ns(t1, t2)) << std::endl;
    std::cout << "peak memory consumption: "
        << format_size(malloc_count_peak() - baseline_memory_alloc) << std::endl;
    std::cout << "compression ratio: " << n / (double) z << std::endl;

    std::cout << std::endl << "running gzip -9" << std::flush;
    std::ofstream file_gz9("fact_gz9");
    malloc_count_reset_peak();
    auto t3 = now();
    lz77::Gzip9Factorizer().factorize(T.begin(), T.end(),
        std::ostream_iterator<lz77::Factor>(file_gz9, ""));
    auto t4 = now();
    uint64_t gz9 = file_gz9.tellp() / sizeof(lz77::Factor);
    file_gz9.close();
    log_runtime(t3, t4);
    std::cout << "throughput: "
        << format_throughput(n, time_diff_ns(t3, t4)) << std::endl;
    std::cout << "peak memory consumption: "
        << format_size(malloc_count_peak() - baseline_memory_alloc) << std::endl;
    std::cout << "compression ratio: " << n / (double) gz9 << std::endl;
}