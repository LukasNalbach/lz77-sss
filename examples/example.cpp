#include <lz77_sss/lz77_sss.hpp>

int main()
{
    // generate a random string
    std::string input = random_repetitive_string(10000, 10000000);
    std::cout << "input size: " << input.size() << std::endl;

    // compute an approximate LZ77 factorization
    auto factorization = lz77_sss<>::factorize_approximate<>(input, {.num_threads = 1});
    std::cout << "number of factors: " << factorization.size() << std::endl;
    std::cout << "compression ratio: " << input.size() / (double)factorization.size() << std::endl;

    // decode the factorization
    auto input_decodeed = lz77_sss<>::decode(factorization, input.size());

    // compute the exact LZ77 factorization to get the approximation ratio
    auto exact_factorization = lz77_sss<>::factorize_exact<>(input, {.num_threads = 1});
    std::cout << "approximation ratio: "
        << factorization.size() / (double)exact_factorization.size() << std::endl;
}