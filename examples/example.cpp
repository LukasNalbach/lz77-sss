#include <lz77_sss/lz77_sss.hpp>

int main()
{
    // generate a random string
    std::string input = random_repetitive_string(10000, 200000);
    std::cout << "input size: " << input.size() << std::endl;
    
    // compute an approximate LZ77 factorization
    std::vector<lz77_sss<>::factor> factorization;
    lz77_sss<>::factorize_approximate<>(input.data(), input.size(), [&](auto f){factorization.emplace_back(f);});
    std::cout << "number of factors: " << factorization.size() << std::endl;
    std::cout << "compression ratio: " << input.size() / (double) factorization.size() << std::endl;

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