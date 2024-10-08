# LZ77-SSS Algorithms
This repository contains implementations of Lempel-Ziv 77 (LZ77) algorithms based on string synchronizing sets [1].

## Usage in C++
### Cmake
```cmake
add_subdirectory(lz77_sss/)
set(LZ77_SSS_BUILD_BENCH_CLI OFF)
set(LZ77_SSS_BUILD_EXAMPLES OFF)
set(LZ77_SSS_BUILD_TESTS OFF)
```

### C++
```c++
#include <lz77_sss/lz77_sss.hpp>
#include <lz77/lpf_factorizer.hpp>

int main() {
    // generate a random string
    std::string input = random_repetitive_string(10000, 10000000);
    std::cout << "input size: " << input.size() << std::endl;

    // compute an approximate LZ77 factorization
    auto factorization = lz77_sss<>::factorize_approximate<>(input);
    std::cout << "number of factors: " << factorization.size() << std::endl;
    std::cout << "compression ratio: " << input.size() / (double) factorization.size() << std::endl;

    // decode the factorization
    auto input_decodeed = lz77_sss<>::decode(factorization, input.size());

    // compute the exact LZ77 factorization to get the approximation ratio
    auto exact_factorization = lz77_sss<>::factorize_exact<>(input);
    std::cout << "approximation ratio: "
        << factorization.size() / (double) exact_factorization.size() << std::endl;
}
```

### References
Jonas Ellert. Sublinear Time Lempel-Ziv (LZ77) Factorization. In String Processing and Information Retrieval (SPIRE) 2023, pages 171-187.