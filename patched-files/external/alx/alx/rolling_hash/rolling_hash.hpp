/*******************************************************************************
 * alx/rolling_hash/rolling_hash.hpp
 *
 * Copyright (C) 2022 Alexander Herlez <alexander.herlez@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#include <assert.h>

#include <bit>
#include <iterator>
#include <random>

#include "rolling_hash/mersenne_modular_arithmetic.hpp"
#include "rolling_hash/modular_arithmetic.hpp"

namespace alx::rolling_hash {
__extension__ typedef unsigned __int128 uint128_t;

template <size_t t_prime_exp = 107>
class rk_prime {
 public:
  rk_prime() : rk_prime(1, 0) {
  }

  rk_prime(uint128_t tau, uint128_t base = 0) : m_tau(tau), m_fp(0) {
    if (base == 0) {
      uint64_t max =
          (t_prime_exp > 64)
              ? ((uint64_t{1} << (127 - std::bit_width(m_prime))) - 1)
              : m_prime - 1;
      m_base = static_cast<uint128_t>(random64(257, max));
    } else {
      m_base = base;
    }
    // prime should be mersenne and (prime*base + prime) should not overflow
    static_assert(t_prime_exp == 107 || t_prime_exp == 61 || t_prime_exp == 89);
    assert(std::bit_width(m_prime) + std::bit_width(m_base) <= 127);

    fill_influence_table();
  }

  // Roll the window by specifying the character that is rolled out of the
  // window and the character that is rolled in the window.
  inline uint128_t roll_in(unsigned char in) {
    return roll(0, in);
  }

  // Roll the window by specifying the character that is rolled out of the
  // window and the character that is rolled in the window.
  inline uint128_t roll_in(uint128_t fp, unsigned char in) const {
    return roll(fp, 0, in);
  }


  // Roll the window by specifying the character that is rolled out of the
  // window and the character that is rolled in the window.
  inline uint128_t roll_out(unsigned char out) {
    return roll(out, 0);
  }


  // Roll the window by specifying the character that is rolled out of the
  // window and the character that is rolled in the window.
  inline uint128_t roll_out(uint128_t fp, unsigned char out) const {
    return roll(fp, out, 0);
  }

  // Roll the window by specifying the character that is rolled out of the
  // window and the character that is rolled in the window.
  inline uint128_t roll(uint128_t fp, unsigned char out, unsigned char in) const {
    fp *= m_base;
    fp = mersenne::mod<uint128_t, m_prime>(fp + m_char_influence[out][in]);
    return fp;
  }

  // Roll the window by specifying the character that is rolled out of the
  // window and the character that is rolled in the window.
  inline uint128_t roll(unsigned char out, unsigned char in) {
    m_fp *= m_base;
    m_fp = mersenne::mod<uint128_t, m_prime>(m_fp + m_char_influence[out][in]);
    return m_fp;
  }

  // Return the prime number used for the rolling hash function.
  inline constexpr uint128_t get_prime() const {
    return m_prime;
  }

  // Return the prime exponent used for the rolling hash function.
  inline constexpr uint128_t get_prime_exp() const {
    return t_prime_exp;
  }

  // Return the fingerprint of the current window.
  inline uint128_t get_fp() const {
    return m_fp;
  }

  void reset() {
    m_fp = 0;
  }

  // Return the base of rolling hash function.
  inline uint128_t get_base() const {
    return m_base;
  }

 private:
  static constexpr uint128_t m_prime = (uint128_t{1} << t_prime_exp) - 1;
  uint128_t m_tau;
  uint128_t m_fp;

  uint128_t m_base;
  uint128_t m_char_influence[256][256] = {};

  // Return a random number that will be used as the base.
  inline static uint64_t random64(uint64_t min, uint64_t max) {
    static std::mt19937_64 g = std::mt19937_64(std::random_device()());
    return (std::uniform_int_distribution<uint64_t>(min, max))(g);
  }

  // Fill up the table needed for fast rolling.
  void fill_influence_table() {
    const uint128_t base_pow_tau_mod_prime =
        modular::pow_mod<uint128_t>(m_base, m_tau, m_prime);
    const uint128_t minus_base_pow_tau_mod_prime =
        mersenne::additive_inverse_mod<uint128_t, m_prime>(
            base_pow_tau_mod_prime);

    // Fill first row
    m_char_influence[0][0] = 0;
    for (size_t j = 1; j < 256; ++j) {
      m_char_influence[0][j] = j;
    }

    // Fill all other rows
    for (size_t i = 1; i < 256; ++i) {
      m_char_influence[i][0] = mersenne::add_mod<uint128_t, m_prime>(
          m_char_influence[i - 1][0], minus_base_pow_tau_mod_prime);
      for (size_t j = 1; j < 256; ++j) {
        m_char_influence[i][j] = mersenne::add_mod<uint128_t, m_prime>(
            m_char_influence[i][j - 1], 1);
      }
    }
  }
};
}  // namespace alx::rolling_hash