# LZ77-SSS-Approximations
This repository contains implementations of Lempel-Ziv 77 (LZ77) approximation algorithms based on string synchronizing sets [1].

## Usage in C++
### Cmake
```cmake
add_subdirectory(lz77_sss_approx/)
set(LZ77_SSS_APRX_BUILD_BENCH_CLI OFF)
set(LZ77_SSS_APRX_BUILD_EXAMPLES OFF)
set(LZ77_SSS_APRX_BUILD_TESTS OFF)
```

### C++
```c++
#include <lz77_sss_approx/lz77_sss_approx.hpp>
#include <lz77/lpf_factorizer.hpp>

int main() {
    // generate a random string
    std::string input = random_repetitive_string(10000000);
    std::cout << "input size: " << input.size() << std::endl;

    // compute the factorization
    auto factorization = lz77_sss_approx<>::factorize(input);
    std::cout << "number of factors: " << factorization.size() << std::endl;
    std::cout << "compression ratio: " << input.size() / (double) factorization.size() << std::endl;

    // decode the factorization
    auto input_decoded = lz77_sss_approx<>::decode(factorization, input.size());

    // compute the exact LZ77 factorization to get the approximation ratio
    std::vector<lz77::Factor> exact_factorization;
    lz77::LPFFactorizer().factorize(input.begin(), input.end(),
        std::back_insert_iterator(exact_factorization));
    std::cout << "approximation ratio: "
        << factorization.size() / (double) exact_factorization.size() << std::endl;
}
```

### References
Jonas Ellert. Sublinear Time Lempel-Ziv (LZ77) Factorization. In String Processing and Information Retrieval (SPIRE) 2023, pages 171-187.