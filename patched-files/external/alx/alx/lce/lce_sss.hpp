/*******************************************************************************
 * alx/lce/lce_sss.hpp
 *
 * Copyright (C) 2022 Alexander Herlez <alexander.herlez@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <chrono>
#include <cmath>
#include <memory>
#include <vector>

#include "lce/lce_classic_for_sss.hpp"
#include "lce/lce_naive_wordwise_xor.hpp"
#include "pred/pred_index.hpp"
#include "rolling_hash/reduce_fingerprints.hpp"
#include "rolling_hash/string_synchronizing_set.hpp"

#ifdef ALX_BENCHMARK_INTERNAL
#include <fmt/core.h>
#include <fmt/ranges.h>

#include "util/timer.hpp"
#ifdef ALX_BENCHMARK_SPACE
#include <malloc_count/malloc_count.h>
#endif
#endif

namespace alx::lce {

template <typename t_char_type = uint8_t, uint64_t t_tau = 1024,
          typename t_index_type = uint32_t, bool t_prefer_long = false>
class lce_sss {
 public:
  typedef t_char_type char_type;
  __extension__ typedef unsigned __int128 uint128_t;
  lce_sss() : m_text(nullptr), m_size(0) {}

  lce_sss(char_type const* text, size_t size) : m_text(text), m_size(size) {
    assert(sizeof(t_char_type) == 1);

#ifdef ALX_BENCHMARK_INTERNAL
    alx::util::timer t;
#ifdef ALX_BENCHMARK_SPACE
    size_t mem_before = malloc_count_current();
    malloc_count_reset_peak();
#endif
#endif

    m_sync_set = rolling_hash::sss<t_index_type, t_tau>(text, size, false);
    // check_string_synchronizing_set(text, m_sync_set);

#ifdef ALX_BENCHMARK_INTERNAL
    fmt::print(" sss_construct_time={}", t.get_and_reset());
    fmt::print(" sss_size={}", m_sync_set.size());
    fmt::print(" sss_runs={}", m_sync_set.num_runs());

#ifdef ALX_BENCHMARK_SPACE
    fmt::print(" sss_construct_mem={}", malloc_count_current() - mem_before);
    fmt::print(" sss_construct_mem_peak={}", malloc_count_peak() - mem_before);
#endif
#endif

#ifdef ALX_BENCHMARK_INTERNAL
#ifdef ALX_BENCHMARK_SPACE
    mem_before = malloc_count_current();
    malloc_count_reset_peak();
#endif
#endif
    m_pred = alx::pred::pred_index<t_index_type, std::bit_width(t_tau) - 1,
                                   t_index_type>(m_sync_set.get_sss());

#ifdef ALX_BENCHMARK_INTERNAL
    fmt::print(" pred_construct_time={}", t.get_and_reset());
#ifdef ALX_BENCHMARK_SPACE
    fmt::print(" pred_construct_mem={}", malloc_count_current() - mem_before);
    fmt::print(" pred_construct_mem_peak={}", malloc_count_peak() - mem_before);
#endif
#endif

#ifdef ALX_BENCHMARK_INTERNAL
#ifdef ALX_BENCHMARK_SPACE
    mem_before = malloc_count_current();
    malloc_count_reset_peak();
#endif
#endif

    std::vector<t_index_type> const& sss = m_sync_set.get_sss();
    std::vector<t_index_type> reduced_fps = reduce_fps_3tau_lexicographic(
        reinterpret_cast<uint8_t const*>(m_text), m_size, m_sync_set);

#ifdef ALX_BENCHMARK_INTERNAL
    fmt::print(" meta_symbols_time={}", t.get_and_reset());
#ifdef ALX_BENCHMARK_SPACE
    fmt::print(" meta_reduction_mem={}", malloc_count_current() - mem_before);
    fmt::print(" meta_reduction_mem_peak={}", malloc_count_peak() - mem_before);
#endif
#endif

#ifdef ALX_BENCHMARK_INTERNAL
#ifdef ALX_BENCHMARK_SPACE
    mem_before = malloc_count_current();
    malloc_count_reset_peak();
#endif
#endif

    m_fp_lce = alx::lce::lce_classic_for_sss<t_index_type, t_tau>(
        reinterpret_cast<uint8_t const*>(m_text), m_size, reduced_fps.data(),
        reduced_fps.size(), sss);

#ifdef ALX_BENCHMARK_INTERNAL
    fmt::print(" meta_lce_construct_time={}", t.get_and_reset());
#ifdef ALX_BENCHMARK_SPACE
    fmt::print(" meta_lce_construct_mem={}",
               malloc_count_current() - mem_before);
    fmt::print(" meta_lce_construct_mem_peak={}",
               malloc_count_peak() - mem_before);
#endif
#endif
  }

  template <typename C>
  lce_sss(C const& container) : lce_sss(container.data(), container.size()) {}

  // Return the number of common letters in text[i..] and text[j..].
  size_t lce(size_t i, size_t j) const {
    if (i == j) [[unlikely]] {
      assert(i < m_size);
      return m_size - i;
    }
    return lce_uneq(i, j);
  }

  // Return the number of common letters in text[i..] and text[j..]. Here i
  // and j must be different.
  size_t lce_uneq(size_t i, size_t j) const {
    assert(i != j);

    size_t l = std::min(i, j);
    size_t r = std::max(i, j);

    return lce_lr(l, r);
  }

  // Return the number of common letters in text[i..] and text[j..].
  // Here l must be smaller than r.
  inline uint64_t lce_lr(size_t l, size_t r) const {
    std::vector<t_index_type> const& sss = m_sync_set.get_sss();
    size_t l_, r_;

    if constexpr (t_prefer_long) {
      // Only scan until synchronizing position
      size_t lce_max{m_size - r};
      size_t lce_local_max{std::min(3 * t_tau, lce_max)};

      pred::result l_res = m_pred.successor(l);
      pred::result r_res = m_pred.successor(r);
      l_ = l_res.pos;
      r_ = r_res.pos;
      if (l_res.exists && r_res.exists && (sss[l_] - l == sss[r_] - r)) {
        lce_local_max =
            std::min(lce_local_max, static_cast<size_t>(sss[l_] - l));
      }

      size_t lce_local = alx::lce::lce_naive_wordwise_xor<t_char_type>::lce_lr(
          m_text, r + lce_local_max, l, r);

      // Case 0: Mismatch at first 3*tau symbols
      if (lce_local < lce_local_max || lce_local == lce_max) {
        return lce_local;
      }
    } else {
      // Naive part until synchronizing position
      size_t lce_max{m_size - r};
      size_t lce_local_max{std::min(3 * t_tau, lce_max)};
      size_t lce_local = alx::lce::lce_naive_wordwise_xor<t_char_type>::lce_lr(
          m_text, r + lce_local_max, l, r);

      // Case 0: Mismatch at first 3*tau symbols
      if (lce_local < lce_local_max || lce_local == lce_max) {
        return lce_local;
      }
      l_ = m_pred.successor(l).pos;
      r_ = m_pred.successor(r).pos;
    }

    if (sss[l_] - l != sss[r_] - r) {
      // Case 1: Positions l' and r' don't sync, (because they are at the end of
      // runs).
      size_t final_lce = std::min(sss[l_] - l, sss[r_] - r) + 2 * t_tau - 1;
      assert(final_lce == alx::lce::lce_naive_wordwise_xor<t_char_type>::lce_lr(
                              m_text, m_size, l, r));
      return final_lce;
    } else {
      // Case 2: Positions l' and r' are synchronized.
      size_t final_lce = (sss[l_] - l) + m_fp_lce.lce_lr(l_, r_);
      assert(final_lce == alx::lce::lce_naive_wordwise_xor<t_char_type>::lce_lr(
                              m_text, m_size, l, r));
      return final_lce;
    }
  }

  // Return {b, lce}, where lce is the number of common letters in text[i..]
  // and text[j..] and b tells whether the lce ends with a mismatch.
  std::pair<bool, size_t> lce_mismatch(size_t i, size_t j) {
    if (i == j) [[unlikely]] {
      assert(i < m_size);
      return {false, m_size - i};
    }

    size_t l = std::min(i, j);
    size_t r = std::max(i, j);

    size_t lce = lce_lr(l, r);
    return {r + lce != m_size, lce};
  }

  // Return whether text[i..] is lexicographic smaller than text[j..]. Here i
  // and j must be different.
  bool is_leq_suffix(size_t i, size_t j) {
    assert(i != j);
    size_t lce_val = lce_uneq(i, j);
    return (
        i + lce_val == m_size ||
        ((j + lce_val != m_size) && m_text[i + lce_val] < m_text[j + lce_val]));
  }

  char_type operator[](size_t i) { return m_text[i]; }

  size_t size() { return m_size; }

  const std::vector<t_index_type>& get_sync_set() {
    return m_sync_set.get_sss();
  }

  const std::vector<uint32_t>& get_ssa() {
    return m_fp_lce.get_ssa();
  }

  const std::vector<uint32_t>& get_issa() {
    return m_fp_lce.get_issa();
  }

  void delete_ssa() {
    return m_fp_lce.delete_ssa();
  }

 private:
  char_type const* m_text;
  size_t m_size;

  alx::pred::pred_index<t_index_type, std::bit_width(t_tau) - 1, t_index_type>
      m_pred;
  rolling_hash::sss<t_index_type, t_tau> m_sync_set;
  alx::lce::lce_classic_for_sss<t_index_type, t_tau> m_fp_lce;
};
}  // namespace alx::lce
/******************************************************************************/
