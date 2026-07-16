# LZ77-SSS Algorithms
This repository contains implementations of Lempel-Ziv 77 (LZ77) algorithms [1] based on string synchronizing sets [2].

## CLI Build Instructions
This implementation has been tested on Ubuntu 24.04 with gcc-12, g++-12, zstd, p7zip-full, gzip, bzip2, xz-utils, lz4, libgtest-dev, libtbb-dev and libomp-dev installed. [bsc](https://github.com/IlyaGrebnov/libbsc) has to be built and installed manually.

```shell
git clone --recurse-submodules https://github.com/LukasNalbach/lz77-sss.git
mkdir build
cd build
cmake ..
cp -rf ../patched-files/* ..
make
```

This creates the following executables in the build/ folder.

## CLI Programs
### Compression Tools
- lz77-sss-3-aprx (LZ77 3-approximation)
- lz77-sss-lpf-lnf-aprx (LZ77 LPF/LNF Approximation)
- lz77-sss-exact (LZ77 Exact Algorithm without interval sampling)
- lz77-sss-exact-smpl (LZ77 Exact Algorithm with interval sampling)
- lz77-sss-decode (decodes/reverts the factorization output by one of the above executables)
- ssszip ((de-)compression tool with an interface that is similar to state-of-the-art compressors)

### Benchmark Tools
- lz77-sss-bench (benchmarks all LZ77 algorithms)
- lz77-sss-bench-tau (benchmarks the LZ77 3-approximation with different values for tau)
- gen-range-queries (generates range query data for a given text)
- bench-range-queries (benchmarks all range data structures using generated query data)
- zip-bench (benchmarks lz4, xz, 7zip, gzip, bzip2, zstd and bsc)

### Examples
- lz77-sss-example (the C++ example from below)

### Test Executables
- test-decomposed-range
- test-dynamic-range
- test-lz77-sss
- test-rabin-karp-substring
- test-sample-index
- test-static-weighted-range

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

int main()
{
    // generate a random string
    std::string input = random_repetitive_string(10000, 200000);
    std::cout << "input length: " << input.size() << std::endl;
    
    // compute an approximate LZ77 factorization
    std::vector<lz77_sss<>::factor> factorization;
    lz77_sss<>::factorize_approximate<>(input.data(), input.size(), [&](auto f){factorization.emplace_back(f);});
    std::cout << "num. of factors: " << factorization.size() << std::endl;
    std::cout << "input length / num. of factors: " << input.size() / (double) factorization.size() << std::endl;

    // decode the factorization
    std::string input_decoded;
    no_init_resize(input_decoded, input.size());
    lz77_sss<>::decode(factorization.begin(), input_decoded.data(), input.size());

    // compute the exact LZ77 factorization to get the approximation ratio
    std::vector<lz77_sss<>::factor> exact_factorization;
    lz77_sss<>::factorize_exact<>(input.data(), input.size(),
        [&](auto f){exact_factorization.emplace_back(f);}, {.num_threads = 1});
    std::cout << "approximation ratio: "
        << factorization.size() / (double) exact_factorization.size() << std::endl;
}
```

## Reproducing the Measurements From the Paper
The measurements in the paper were performed on the following three 50 GiB inputs:

| File | Description | Link |
|---|---|---|
| chr19.50Gi | repeated human chromosome 19 | [drive.google.com](https://drive.google.com/file/d/1kErl5gp-e509Y34tLbrsRB4wDfM2gPHu/view?usp=sharing) |
| sars2.50Gi | SARS-CoV-2 genomes | [drive.google.com](https://drive.google.com/file/d/132omEsl7xcjNuHaWiDd6HM233SWzwFcU/view?usp=sharing) |
| dewiki.50Gi | German Wikipedia dump | [drive.google.com](https://drive.google.com/file/d/1y1ajO7JI3QVQAg_ud1rXkgRgihy0laQl/view?usp=sharing) |

Download the three files (e.g. into a `texts/` folder), for example using [gdown](https://github.com/wkentaro/gdown):
```shell
pip install gdown
mkdir -p texts
gdown 1kErl5gp-e509Y34tLbrsRB4wDfM2gPHu -O texts/chr19.50Gi
gdown 132omEsl7xcjNuHaWiDd6HM233SWzwFcU -O texts/sars2.50Gi
gdown 1y1ajO7JI3QVQAg_ud1rXkgRgihy0laQl -O texts/dewiki.50Gi
```

Then build the CLI programs as described in [CLI Build Instructions](#cli-build-instructions) and run the benchmarks with the `lz77-sss-bench` and `zip-bench` executables in the `build/bench/` folder.

### lz77-sss-bench
Benchmarks all LZ77 algorithms on a given text.
```shell
build/bench/lz77-sss-bench <file> <max_threads> <result_file>
```
- `<file>`: input text (e.g. `texts/chr19.50Gi`)
- `<max_threads>` (optional): maximum number of threads to use (default: all available)
- `<result_file>` (optional): file the results are appended to

### zip-bench
Benchmarks lz4, xz, 7zip, gzip, bzip2, zstd and bsc on a given text.
```shell
build/bench/zip-bench <input_file> <min_threads> <max_threads> <result_file>
```
- `<input_file>`: input text (e.g. `texts/chr19.50Gi`)
- `<min_threads>`: minimum number of threads to use
- `<max_threads>`: maximum number of threads to use
- `<result_file>` (optional): file the results are appended to

### Running all measurements
To reproduce the measurements for all three texts using all available threads and collecting the results in `results.txt`:
```shell
for text in texts/chr19.50Gi texts/sars2.50Gi texts/dewiki.50Gi; do
    build/bench/lz77-sss-bench "$text" $(nproc) results.txt
    build/bench/zip-bench "$text" 1 $(nproc) results.txt
done
```

### References
[1] Jonas Ellert. Sublinear Time Lempel-Ziv (LZ77) Factorization. In String Processing and Information Retrieval (SPIRE) 2023, pages 171-187. ([springer.com](https://link.springer.com/chapter/10.1007/978-3-031-43980-3_14))

[2] Dominik Kempa and Tomasz Kociumaka. String synchronizing sets: sublinear-time BWT construction and optimal LCE data structure. In Proceedings of the 51st Annual ACM SIGACT Symposium on Theory of Computing (STOC) 2019, pages 756-767. ([arxiv.org](https://arxiv.org/abs/1904.04228))