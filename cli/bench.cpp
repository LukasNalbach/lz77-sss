#include <fstream>

#include <lz77_sss_approx/lz77_sss_approx.hpp>
#include <lz77/lpf_factorizer.hpp>
#include <lz77/gzip9_factorizer.hpp>

template <typename pos_t>
void check_correctness(std::ifstream& file, std::string& T) {
    std::string T_rev = lz77_sss_approx<pos_t>::decode(file,T.size());
    bool correct = true;

    if (T_rev.size() != T.size()) {
        correct = false;
    } else {
        for (pos_t i=0; i<T.size(); i++) {
            if (T_rev[i] != T[i]) {
                correct = false;
                break;
            }
        }
    }

    std::cout << "the factorization is " << (correct ? "" : "not ") << "correct" << std::endl;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "usage: bench <file>";
        exit(-1);
    }

    std::ifstream input_file;
    input_file.open(argv[1]);
    
    if (!input_file.good()) {
        std::cout << "error: could not read <file>";
        exit(-1);
    }

    input_file.seekg(0,std::ios::end);
    uint64_t n = input_file.tellg();
    input_file.seekg(0,std::ios::beg);
    auto time = now();
    std::cout << "reading T (" << format_size(n) << ")" << std::flush;
    std::string T;
    no_init_resize(T,n);
    read_from_file(input_file,(const char*)&T[0],n);
    input_file.close();
    time = log_runtime(time);
    std::cout << std::endl << "running LZ77 approximation" << std::endl;

    {
        std::ofstream file_aprx("fact_aprx");

        if (n <= UINT_MAX) {
            lz77_sss_approx<uint32_t>::factorize<greedy,lpf_optimal>(T,file_aprx,true);
            file_aprx.close();
            std::ifstream fact_aprx_file("fact_aprx");
            check_correctness<uint32_t>(fact_aprx_file,T);
        } else {
            lz77_sss_approx<uint64_t>::factorize<greedy,lpf_optimal>(T,file_aprx,true);
            file_aprx.close();
            std::ifstream fact_aprx_file("fact_aprx");
            check_correctness<uint64_t>(fact_aprx_file,T);
        }
    }

    std::cout << std::endl;
    std::cout << "computing LZ77 using the LPF array" << std::flush;
    uint64_t baseline_memory_alloc = malloc_count_current();

    {
        std::ofstream file_exact("fact_exact");
        malloc_count_reset_peak();

        time = now();
        lz77::LPFFactorizer().factorize(T.begin(),T.end(),std::ostream_iterator<lz77::Factor>(file_exact,""));
        time = log_runtime(time);

        uint64_t z = file_exact.tellp()/sizeof(lz77::Factor);
        file_exact.close();

        std::cout << "peak memory consumption: " << format_size(malloc_count_peak()-baseline_memory_alloc) << std::endl;
        std::cout << "compression ratio: " << n/(double)z << std::endl;
    }

    std::cout << std::endl;
    std::cout << "computing gzip -9 output" << std::flush;

    {
        std::ofstream file_gz9("fact_gz9");

        malloc_count_reset_peak();
        time = now();
        lz77::Gzip9Factorizer().factorize(T.begin(),T.end(),std::ostream_iterator<lz77::Factor>(file_gz9,""));
        time = log_runtime(time);

        uint64_t gz9 = file_gz9.tellp()/sizeof(lz77::Factor);
        file_gz9.close();
        
        std::cout << "peak memory consumption: " << format_size(malloc_count_peak()-baseline_memory_alloc) << std::endl;
        std::cout << "compression ratio: " << n/(double)gz9 << std::endl;
    }
}