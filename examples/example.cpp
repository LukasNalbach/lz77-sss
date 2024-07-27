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
    auto input_decodeed = lz77_sss_approx<>::decode(factorization, input.size());

    // compute the exact LZ77 factorization to get the approximation ratio
    std::vector<lz77::Factor> exact_factorization;
    lz77::LPFFactorizer().factorize(input.begin(), input.end(),
        std::back_insert_iterator(exact_factorization));
    std::cout << "approximation ratio: "
        << factorization.size() / (double) exact_factorization.size() << std::endl;
}