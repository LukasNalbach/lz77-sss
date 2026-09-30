# LZ77-SSS Algorithms
This repository contains implementations of Lempel-Ziv 77 (LZ77) algorithms [1] based on string synchronizing sets [2]. These implementations are described in [3] (accepted at ALENEX 2027, [arxiv.org](https://arxiv.org/abs/2609.30193)).

## CLI Build Instructions
This implementation has been tested on Ubuntu 26.04 with g++ 15.2, CMake 4.2, git and libtbb-dev installed. ssszip needs the chosen postcompressor in the PATH (by default [bsc](https://github.com/IlyaGrebnov/libbsc), which has to be built and installed manually); zip-bench needs time and taskset and benchmarks [alz](https://github.com/pdinklag/alz) (built automatically from external/alz, followed by bsc) and those of lz4, xz, 7z, gzip, bzip2, zstd and bsc that are installed.

```shell
git clone --recurse-submodules https://github.com/LukasNalbach/lz77-sss.git
mkdir build
cd build
cmake ..
make
```

cmake applies the patches in patches/ to the lce submodule with git apply.

### CMake Options
| Option | Default | Description |
|---|---|---|
| LZ77_SSS_BUILD_CLI | ON | build the CLI programs |
| LZ77_SSS_BUILD_BENCH | ON (OFF on Windows) | build the benchmark programs and alz |
| LZ77_SSS_BUILD_EXAMPLES | OFF | build the example program |
| LZ77_SSS_BUILD_TESTS | ON | build the tests |
| LZ77_SSS_USE_MALLOC_COUNT | ON | track memory usage via malloc_count (overrides malloc/free) |
| LZ77_SSS_MARCH_NATIVE | ON | tune the build for the host CPU (-march=native); turn OFF for portable binaries |
| LZ77_SSS_ENABLE_LTO | ON | enable link-time optimization for the executables |
| LZ77_SSS_EXTERN_TEMPLATES | ON | compile the factorizer once for the plain, packed and split text instead of in every translation unit |
| LZ77_SSS_LCE_DIR | external/lce | lce source directory |
| LZ77_SSS_WINDOWS_TBB_ROOT | | Windows only: root of an extracted oneAPI TBB release |

This creates the compression tools in `build/cli/`, the benchmark tools in `build/bench/`,
the examples in `build/examples/` and the tests in `build/tests/`.

## CLI Programs
### Compression Tools
- lz77-sss-3-aprx (LZ77 3-approximation)
- lz77-sss-exact (LZ77 Exact Algorithm without interval sampling)
- lz77-sss-exact-smpl (LZ77 Exact Algorithm with interval sampling)
- lz77-sss-decode (decodes/reverts the factorization output by one of the above executables)
- ssszip ((de-)compression tool with an interface that is similar to state-of-the-art compressors)

### Benchmark Tools
- lz77-sss-bench (benchmarks all LZ77 algorithms)
- lz77-sss-bench-tau (benchmarks the LZ77 3-approximation with different values for tau)
- gen-range-queries (generates range query data for a given text)
- bench-range-queries (benchmarks all range data structures using generated query data)
- zip-bench (benchmarks alz, lz4, xz, 7zip, gzip, bzip2, zstd and bsc)

### Examples
- lz77-sss-example (the C++ example from below, built with LZ77_SSS_BUILD_EXAMPLES=ON)

### Test Executables
- test-decomposed-range
- test-dynamic-range
- test-lz77-sss
- test-rabin-karp-substring
- test-sample-index
- test-static-weighted-range

## Compression Tool Interface
```
usage: lz77-sss-3-aprx [options] <input_file>
 -o <output_file>  output file (default: <input_file>.lz77-aprx)
 -t <threads>      number of threads (default: all)
 -enc <encoding>   how the text is kept in memory (default: auto):
                   plain   one byte per character (fastest)
                   packed  fewer bits per character for small alphabets
                   split   fewer bits for frequent characters (slowest)
                   auto    packed or split if that saves memory, else plain
 -fasta <mode>     handle the header lines of FASTA files separately from
                   the sequences: on, off or auto (default: auto, on for
                   FASTA files)
 -h                show help
```
The same interface is used by lz77-sss-exact and lz77-sss-exact-smpl; their default output
file is <input_file>.lz77-exact. The exact algorithms additionally accept the following options:
```
 -alg <algorithm>  exact algorithm (default: auto):
                   sss     string synchronizing sets, for repetitive texts
                   sa      suffix array, needs about 9 bytes per character
                   auto    sa if it needs at most 25 % more memory than sss,
                           else sss
 -r <range_ds>     range data structure (default: d-swsg):
                   swsg    static weighted square grid
                   swkdt   static weighted k-d tree
                   sdsg    semi-dynamic square grid (single-threaded)
                   d-...   one structure per character, e.g. d-swsg
```

The input text is scanned in parallel before the factorization to compute its character
histogram, from which the default encoding is chosen as described for -enc.

## Decoding Tool Interface
```
usage: lz77-sss-decode [options] <input_file> <output_file>
 -ram          keep the whole output in memory while decoding (default: only
               the last 64 MiB, older parts are read back from <output_file>)
 -t <threads>  number of threads (default: all)
 -h            show help
```

## ssszip Interface
```
usage: ssszip [options] <input_file>
 -d                    decompress <input_file> (*.ssszip.<postcompressor>)
 -ram                  keep the whole output in memory while decompressing
                       (default: only the last 64 MiB)
 -o <base_name>        write <base_name>.ssszip.<postcompressor> (default:
                       <input_file>); with -d: the output file (default:
                       <input_file> without .ssszip.<postcompressor>)
 -t <threads>          number of threads (default: all)
 -pc <postcompressor>  postcompressor: bsc (default), zstd, xz, lzma, gzip,
                       pigz, bzip2, pbzip2, lbzip2, brotli, lz4, lzop, lzip,
                       plzip, bzip3, 7z
 -<quality>            post-compression quality, e.g. -9: higher gives smaller
                       output but takes longer (range and default depend on
                       the postcompressor); for bsc and bzip3 the block size
                       in MB (bsc default: largest block that needs no more
                       memory than the factorization, at most 2047)
 -enc <encoding>       how the text is kept in memory (default: auto):
                       plain   one byte per character (fastest)
                       packed  fewer bits per character for small alphabets
                       split   fewer bits for frequent characters (slowest)
                       auto    packed or split if that saves memory, else plain
 -fasta <mode>         store the header lines of FASTA files separately from
                       the sequences: on, off or auto (default: auto, on for
                       FASTA files)
 -q                    print nothing
 -v                    print progress and statistics (default)
 -m <m_file>           append results to <m_file>
 -h                    show help
```

## Benchmark Tool Interfaces
```
usage: lz77-sss-bench [-m <m_file>] <input_file> [<max_threads>]
 -m <m_file>    append results to <m_file>
 <max_threads>  largest thread count; the count doubles from 1 (default: all)

usage: lz77-sss-bench-tau [-m <m_file>] <input_file> [<min_tau> <max_tau>]
 -m <m_file>  append results to <m_file>
 <min_tau>    smallest tau (default: 4, at least 4)
 <max_tau>    largest tau (default: 4096, at most 4096); every power of two
              in between is run

usage: gen-range-queries <input_file> <queries_file>

usage: bench-range-queries [-m <m_file>] <input_file> <queries_file>
                           [<min_win> <max_win>]
 -m <m_file>     append results to <m_file>
 <queries_file>  queries written by gen-range-queries for <input_file>
 <min_win>       log2 of the smallest square grid window (default: 11)
 <max_win>       log2 of the largest square grid window (default: 16)

usage: zip-bench [-m <m_file>] <input_file> <min_threads> <max_threads>
 -m <m_file>    append results to <m_file>
 <min_threads>  smallest thread count for the multithreaded compressors
 <max_threads>  largest thread count; the count doubles from <min_threads>
```
With -m, each benchmark appends one `RESULT key=value ...` line per run to <m_file>.

## Usage in C++
### Cmake
```cmake
set(LZ77_SSS_BUILD_CLI OFF CACHE BOOL "" FORCE)
set(LZ77_SSS_BUILD_BENCH OFF CACHE BOOL "" FORCE)
set(LZ77_SSS_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(LZ77_SSS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory(lz77_sss/)
target_link_libraries(<your_target> PRIVATE lz77_sss)
```

### C++
```c++
#include <lz77_sss/inst.hpp>
#include <lz77_sss/misc/repetitive_input.hpp>

int main()
{
    // generate a random string
    std::string input = random_repetitive_input<std::string>(4'000, 1'000'000);
    std::cout << "input length: " << input.size() << std::endl;
    
    // compute an approximate LZ77 factorization
    std::vector<lz77_sss::factor> factorization;
    lz77_sss::factorize_approximate<>(input.data(), input.size(), [&](auto f){factorization.emplace_back(f);});
    std::cout << "num. of factors: " << factorization.size() << std::endl;
    std::cout << "input length / num. of factors: " << input.size() / (double) factorization.size() << std::endl;

    // decode the factorization
    std::string input_decoded;
    no_init_resize(input_decoded, input.size());
    lz77_sss::decode(factorization.begin(), input_decoded.data(), input.size());

    // compute the exact LZ77 factorization to get the approximation ratio
    std::vector<lz77_sss::factor> exact_factorization;
    lz77_sss::factorize_exact<>(input.data(), input.size(),
        [&](auto f){exact_factorization.emplace_back(f);});
    std::cout << "approximation ratio: "
        << factorization.size() / (double) exact_factorization.size() << std::endl;
}
```

### Integer Alphabets
`factorize_approximate` and `factorize_exact` also accept a `const uint32_t*` input. They replace
each value by its rank in the sorted alphabet, keep the ranks bit-packed and report literal factors
with the original values. A text that already consists of ranks in `[0, sigma)` can be passed as one
of these text types:

| type                        | space per symbol                  |
|-----------------------------|-----------------------------------|
| `lz77_sss::int_direct_text` | 4 bytes (the text is not copied)  |
| `lz77_sss::int_packed_text` | ⌈log2 sigma⌉ bits                 |
| `lz77_sss::int_split_text`  | fewer bits for frequent symbols   |

```c++
std::vector<uint32_t> ranks = ...; // all values < sigma
std::vector<lz77_sss::factor> factorization;
lz77_sss::factorize_exact(lz77_sss::int_packed_text(ranks.data(), ranks.size(), sigma),
    [&](auto f){factorization.emplace_back(f);});
std::vector<uint32_t> decoded(ranks.size());
lz77_sss::decode(factorization.begin(), decoded.data(), decoded.size());
```

## References
[1] Jonas Ellert. Sublinear Time Lempel-Ziv (LZ77) Factorization. In String Processing and Information Retrieval (SPIRE) 2023, pages 171-187. ([springer.com](https://link.springer.com/chapter/10.1007/978-3-031-43980-3_14))

[2] Dominik Kempa and Tomasz Kociumaka. String synchronizing sets: sublinear-time BWT construction and optimal LCE data structure. In Proceedings of the 51st Annual ACM SIGACT Symposium on Theory of Computing (STOC) 2019, pages 756-767. ([arxiv.org](https://arxiv.org/abs/1904.04228))

[3] Jonas Ellert and Lukas Nalbach. Practical and Space-Efficient LZ77 and LZ Pre-Compression via String Synchronizing Sets. Accepted at the SIAM Symposium on Algorithm Engineering and Experiments (ALENEX) 2027. ([arxiv.org](https://arxiv.org/abs/2609.30193))
