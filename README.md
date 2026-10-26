# LZ77-SSS

This repository contains parallel algorithms that compute the Lempel-Ziv 77 (LZ77) factorization [1] of a text with
string synchronizing sets [2], described in [3], together with command-line tools, benchmarks and tests.

- `lz77_sss::factorize_approximate` computes a 3-approximation of the LZ77 factorization (at most 3 times as many
  factors), `lz77_sss::factorize_exact` the exact factorization and `lz77_sss::decode` the text of a factorization.
- The text consists of bytes or of integers of 1, 4 or 8 bytes; it can be stored with fewer bits per symbol with the
  text classes of [compact-texts](https://github.com/LukasNalbach/compact-texts).
- **ssszip** compresses files with an LZ77 factorization and a general-purpose postcompressor.

## Contents

- [Build](#build)
- [Command-line tools](#command-line-tools)
- [C++ usage](#c-usage)
- [Tests](#tests)
- [References](#references)
- [License](#license)

## Build

Requirements: a C++23 compiler, CMake 3.21 or later, git and Intel TBB (`libtbb-dev`). Tested on Ubuntu 26.04 with
GCC 15.2 and CMake 4.2.

```shell
git clone https://github.com/LukasNalbach/lz77-sss.git
cd lz77-sss
mkdir build
cd build
cmake ..
make
```

The command-line tools are written to `build/cli/`, the benchmark tools to `build/bench/` and the tests to
`build/tests/`. `cmake` downloads the submodules [lce](https://github.com/LukasNalbach/lce) (with the patch in
[patches/](patches/)), [compact-texts](https://github.com/LukasNalbach/compact-texts),
[string-utils](https://github.com/LukasNalbach/string-utils) (build settings, logging, memory helpers, ips4o, libsais,
malloc_count and googletest) and, for the benchmarks, [alz](https://github.com/pdinklag/alz).

| Option | Description |
| --- | --- |
| `LZ77_SSS_BUILD_CLI` | build the command-line tools (`ON`) |
| `LZ77_SSS_BUILD_BENCH` | build the benchmark tools and alz (`OFF`; not on Windows) |
| `LZ77_SSS_BUILD_TESTS` | build the tests (`ON`) |
| `LZ77_SSS_ENABLE_LTO` | link-time optimization of the executables (`ON`) |
| `LZ77_SSS_LCE_DIR` | the directory of the lce sources (default: `external/lce`) |

The other options (`LZ77_SSS_MARCH_NATIVE`, `LZ77_SSS_USE_MALLOC_COUNT`, ...) and the environment variables are
described in string-utils (with the prefix `LZ77_SSS`).

## Command-line tools

The compression tools print their options with `-h` and take the number of threads with `-t`.

- **lz77-sss-3-aprx**, **lz77-sss-exact**, **lz77-sss-exact-smpl**: write the 3-approximation or the exact
  factorization of a file (the latter without or with interval sampling) to `<input_file>.lz77-aprx` or
  `<input_file>.lz77-exact`. `-enc` selects how the text is stored in memory and `-fasta` stores the header lines of a
  FASTA file separately. `lz77-sss-exact -t 8 text.txt`
- **lz77-sss-decode**: writes the text of a factorization. `lz77-sss-decode text.txt.lz77-exact text.out`
- **ssszip**: compresses a file into `<input_file>.ssszip.<postcompressor>` and decompresses it with `-d`. `-pc` selects
  the postcompressor (default: [bsc](https://github.com/IlyaGrebnov/libbsc)), which has to be installed.
  `ssszip -t 8 text.txt`
- **lz77-sss-bench**, **lz77-sss-bench-tau**, **gen-range-queries**, **bench-range-queries**, **zip-bench** (in
  `build/bench/`): measure the LZ77 algorithms, the 3-approximation for different tau, the range data structures, and
  alz, lz4, xz, 7z, gzip, bzip2, zstd and bsc. With `-m <file>`, they append `RESULT key=value ...` lines to a file.
  `lz77-sss-bench -m results.txt text.txt 16`

## C++ usage

Add lz77-sss to a CMake project and link the target `lz77_sss`:

```cmake
set(CMAKE_CXX_STANDARD 23)
set(LZ77_SSS_BUILD_CLI OFF CACHE BOOL "" FORCE)
set(LZ77_SSS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory(lz77_sss/ EXCLUDE_FROM_ALL)
target_link_libraries(<your_target> PRIVATE lz77_sss)
```

```c++
#include <lz77_sss/inst.hpp>
#include <misc/repetitive_input.hpp>

int main()
{
    // generate a repetitive string
    std::mt19937_64 random(std::random_device {}());
    std::string input = generate_repetitive_input<std::string>(
        random_repetitive_params(1'000'000, random));
    // the algorithms read up to 4 tau bytes after the input
    resize_with_zero_excess(input, input.size(), 4 * lz77_sss::default_tau);
    std::cout << "input length: " << input.size() << std::endl;

    // compute an approximate LZ77 factorization
    std::vector<lz77_sss::factor> factorization;
    lz77_sss::factorize_approximate<>(input.data(), input.size(),
        [&](auto f) { factorization.emplace_back(f); });
    std::cout << "num. of factors: " << factorization.size() << std::endl;
    std::cout << "input length / num. of factors: "
        << input.size() / (double) factorization.size() << std::endl;

    // decode the factorization
    std::string input_decoded(input.size(), '\0');
    lz77_sss::decode(factorization.begin(), input_decoded.data(),
        input.size());
    std::cout << "decoded correctly: " << (input_decoded == input)
        << std::endl;

    // compute the exact LZ77 factorization to get the approximation ratio
    std::vector<lz77_sss::factor> exact_factorization;
    lz77_sss::factorize_exact<>(input.data(), input.size(),
        [&](auto f) { exact_factorization.emplace_back(f); });
    std::cout << "approximation ratio: "
        << factorization.size() / (double) exact_factorization.size()
        << std::endl;
}
```

`factorize_approximate` and `factorize_exact` also take a pointer to signed or unsigned integers of 1, 4 or 8 bytes, or
one of the text types `lz77_sss::int_direct_text` (4 bytes per symbol), `lz77_sss::int64_direct_text` (8 bytes),
`lz77_sss::int_packed_text` (ceil(log2 sigma) bits) and `lz77_sss::int_split_text` (fewer bits for frequent symbols).
If `sigma` is given, all values have to be smaller than `sigma`. The exact factorization may change the input while it
runs and restores it before it returns.

```c++
#include <lz77_sss/inst.hpp>

int main()
{
    // a text of 32-bit integers
    std::mt19937_64 random(std::random_device {}());
    std::vector<uint32_t> input(1'000'000);
    for (uint32_t& c : input) c = random() % 1'000;

    std::vector<lz77_sss::factor> factorization;
    lz77_sss::factorize_exact(input.data(), input.size(),
        [&](auto f) { factorization.emplace_back(f); });
    std::cout << "num. of factors: " << factorization.size() << std::endl;

    std::vector<uint32_t> decoded(input.size());
    lz77_sss::decode(factorization.begin(), decoded.data(), decoded.size());
    std::cout << "decoded correctly: " << (decoded == input) << std::endl;
}
```

## Tests

```shell
build/tests/test-lz77-sss
```

The other tests in `build/tests/` check the range data structures, the Rabin-Karp fingerprints and the sample index.
They check random inputs: at most `LZ77_SSS_TEST_SECONDS` seconds (default 10) or `LZ77_SSS_TEST_ITERATIONS` inputs
(default 1000) per function. A failed test prints its seed; `LZ77_SSS_TEST_SEED=<seed>` repeats its inputs.

## References

[1] Jonas Ellert. Sublinear Time Lempel-Ziv (LZ77) Factorization. In String Processing and Information Retrieval
(SPIRE) 2023, pages 171-187. ([springer.com](https://link.springer.com/chapter/10.1007/978-3-031-43980-3_14))

[2] Dominik Kempa and Tomasz Kociumaka. String synchronizing sets: sublinear-time BWT construction and optimal LCE data
structure. In Proceedings of the 51st Annual ACM SIGACT Symposium on Theory of Computing (STOC) 2019, pages 756-767.
([arxiv.org](https://arxiv.org/abs/1904.04228))

[3] Jonas Ellert and Lukas Nalbach. Practical and Space-Efficient LZ77 and LZ Pre-Compression via String Synchronizing
Sets. Accepted at the SIAM Symposium on Algorithm Engineering and Experiments (ALENEX) 2027.
([arxiv.org](https://arxiv.org/abs/2609.30193))

## License

This project is licensed under the MIT License, see [LICENSE](LICENSE). The submodules in [external/](external/) have
their own licenses.
