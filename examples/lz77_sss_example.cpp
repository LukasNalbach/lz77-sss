/**
 * part of LukasNalbach/lz77-sss
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <lz77_sss/lz77_sss.hpp>

int main()
{
    // generate a random string
    std::string input = random_repetitive_string(4'000, 1'000'000);
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
        [&](auto f){exact_factorization.emplace_back(f);});
    std::cout << "approximation ratio: "
        << factorization.size() / (double) exact_factorization.size() << std::endl;
}