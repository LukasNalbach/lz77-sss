#include <lz77_sss/lz77_sss.hpp>

int main()
{
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