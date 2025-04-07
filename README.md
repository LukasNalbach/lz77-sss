# LZ77-SSS Algorithms
This repository contains implementations of Lempel-Ziv 77 (LZ77) algorithms [1] based on string synchronizing sets [2].

## CLI Build Instructions
This implementation has been tested on Ubuntu 24.04 with gcc-12, g++-12, zstd, p7zip-full, gzip, bzip2, xz-utils, lz4, libtbb-dev and libomp-dev installed. [bsc](https://github.com/IlyaGrebnov/libbsc) has to be built and installed manually.

```shell
clone https://github.com/LukasNalbach/lz77-sss.git
mkdir build
cd build
cmake ..
cp -rf ../patched-files/* ..
make
```

This creates the following executables in the build/ folder.

## CLI Programs
### Compression Tools
- lz77_sss_3-aprx (LZ77 3-Approximation)
- lz77_sss_1_5-aprx (LZ77 LPF/LNF Approximation)
- lz77_sss_exact (LZ77 Exact Algorithm without interval sampling)
- lz77_sss_exact-smpl (LZ77 Exact Algorithm with interval sampling)
- lz77_sss_decode (decodes/reverts the factorization output by one of the above executables)
- ssszip ((de-)compression tool with an interface that is similar to state-of-the-art compressors)

### Benchmark Tools
- lz77_sss_bench (benchmarks all LZ77 algorithms)
- lz77_sss_bench-tau (benchmarks the LZ77 3-Approximation with different values for tau)
- gen_range_queries (generates range query data for a given text)
- bench_range_queries (benchmarks all range data structures using generated query data)
- zip-bench (benchmarks lz4, xz, 7zip, gzip, bzip2, zstd and bsc)

### Examples
- lz77_sss_example (the C++ example from below)

### Test Executables
- test-decomposed_range
- test-dynamic_range
- test-lz77_sss
- test-rk61_substring
- test-sample_index
- test-static_weighted_range

## ssszip Interface
```
usage: ssszip [...] <input_file>
 -d                decompress <input_file> (with extension .ssszip.<encoder>)
 -o <output_file>  output file path (default: <input_file>.ssszip.<encoder>)
 -t <threads>      number of threads to use (default: all)
 -e <encoder>      name of the encoder binary (default: zstd)
 -0/-1/-2/...      encoding quality (default: 4)
 -k                keep (don't delete) <input file>
 -q                quiet mode (disables all logs)
 -v                shows verbose information
 -r <result_file>  write results to <result_file>
 -h                show help
```

## Usage in C++
### Cmake
```cmake
add_subdirectory(lz77_sss/)
set(LZ77_SSS_BUILD_CLI OFF)
set(LZ77_SSS_BUILD_BENCH_CLI OFF)
set(LZ77_SSS_BUILD_EXAMPLES OFF)
set(LZ77_SSS_BUILD_TESTS OFF)
```

### C++
```c++
#include <lz77_sss/lz77_sss.hpp>

int main() {
    // generate a random string
    char* input = nullptr;
    uint32_t input_size = random_repetitive_string(input, 10000, 200000);
    std::cout << "input size: " << input_size << std::endl;
    
    // compute an approximate LZ77 factorization
    std::vector<lz77_sss<>::factor> factorization;
    lz77_sss<>::factorize_approximate<>(input, input_size, [&](auto f){factorization.emplace_back(f);});
    std::cout << "number of factors: " << factorization.size() << std::endl;
    std::cout << "compression ratio: " << input_size / (double) factorization.size() << std::endl;

    // decode the factorization
    char* input_decoded = (char*) std::aligned_alloc(16, input_size);
    lz77_sss<>::decode(factorization.begin(), input_decoded, input_size);

    // compute the exact LZ77 factorization to get the approximation ratio
    std::vector<lz77_sss<>::factor> exact_factorization;
    lz77_sss<>::factorize_exact<>(input, input_size,
        [&](auto f){exact_factorization.emplace_back(f);}, {.num_threads = 1});
    std::cout << "approximation ratio: "
        << factorization.size() / (double) exact_factorization.size() << std::endl;
}
```

### References
[1] Jonas Ellert. Sublinear Time Lempel-Ziv (LZ77) Factorization. In String Processing and Information Retrieval (SPIRE) 2023, pages 171-187. ([springer.com](https://link.springer.com/chapter/10.1007/978-3-031-43980-3_14))

[2] Dominik Kempa and Tomasz Kociumaka. String synchronizing sets: sublinear-time BWT construction and optimal LCE data structure. In Proceedings of the 51st Annual ACM SIGACT Symposium on Theory of Computing (STOC) 2019, pages 756-767. ([arxiv.org](https://arxiv.org/abs/1904.04228))