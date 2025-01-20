/*******************************************************************************
 * lce/ds/lce_classic.hpp
 *
 * Copyright (C) 2022 Alexander Herlez <alexander.herlez@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#include <assert.h>

#include <cstdint>
#include <gsaca-double-sort-par.hpp>

#include "ds/lce_naive_wordwise_xor.hpp"
#include "rmq/rmq_n.hpp"

#ifdef LCE_BENCHMARK_INTERNAL
#include <fmt/core.h>
#include <fmt/ranges.h>

#include "util/timer.hpp"
#ifdef LCE_BENCHMARK_SPACE
#include <malloc_count/malloc_count.h>
#endif
#endif

namespace lce::ds {

template <typename t_char_type = uint8_t, typename t_index_type = uint32_t>
class lce_classic {
 public:
  typedef t_char_type char_type;

  lce_classic() : m_text(nullptr), m_size(0) {
  }

  lce_classic(char_type const* text, size_t size) : m_text(text), m_size(size) {
    std::vector<t_index_type> sa(size);
    // sort sa
#ifdef LCE_BENCHMARK_INTERNAL
    lce::util::timer t;
#ifdef LCE_BENCHMARK_SPACE
    size_t mem_before = malloc_count_current();
    malloc_count_reset_peak();
#endif
#endif
    gsaca_for_lce(text, sa.data(), size);
#ifdef LCE_BENCHMARK_INTERNAL
    fmt::print(" sa_time={}", t.get_and_reset());
#ifdef LCE_BENCHMARK_SPACE
    fmt::print(" sa_mem={}", malloc_count_current() - mem_before);
    fmt::print(" sa_mem_peak={}", malloc_count_peak() - mem_before);
    mem_before = malloc_count_current();
    malloc_count_reset_peak();
#endif
#endif

    // build isa
    m_isa.resize(size);
#pragma omp parallel for
    for (size_t i = 0; i < sa.size(); ++i) {
      m_isa[sa[i]] = i;
    }

#ifdef LCE_BENCHMARK_INTERNAL
    fmt::print(" isa_time={}", t.get_and_reset());
#ifdef LCE_BENCHMARK_SPACE
    fmt::print(" isa_mem={}", malloc_count_current() - mem_before);
    fmt::print(" isa_mem_peak={}", malloc_count_peak() - mem_before);
    mem_before = malloc_count_current();
    malloc_count_reset_peak();
#endif
#endif

    // build lcp
    m_lcp.resize(sa.size());
    m_lcp[0] = 0;
    size_t current_lcp = 0;

#pragma omp parallel
    {
      const int t = omp_get_thread_num();
      const int nt = omp_get_num_threads();
      const size_t slice_size = size / nt;

      const size_t begin = t * slice_size;
      const size_t end = (t < nt - 1) ? (t + 1) * slice_size : size;

      size_t current_lcp = 0;
      for (size_t i{begin}; i < end; ++i) {
        size_t suffix_array_pos = m_isa[i];
        if (suffix_array_pos == 0) {
          continue;
        }
        assert(suffix_array_pos != 0);

        size_t preceding_suffix_pos = sa[suffix_array_pos - 1];
        current_lcp += lce_naive_wordwise_xor<char_type>::lce_uneq(
            text, size, i + current_lcp, preceding_suffix_pos + current_lcp);
        m_lcp[suffix_array_pos] = current_lcp;
        assert(lce_naive_wordwise_xor<char_type>::lce_uneq(
                   text, size, i, preceding_suffix_pos) == current_lcp);

        if (current_lcp != 0) {
          --current_lcp;
        }
      }
    }

#ifdef LCE_BENCHMARK_INTERNAL
    fmt::print(" lcp_time={}", t.get_and_reset());
#ifdef LCE_BENCHMARK_SPACE
    fmt::print(" lcp_mem={}", malloc_count_current() - mem_before);
    fmt::print(" lcp_mem_peak={}", malloc_count_peak() - mem_before);
    mem_before = malloc_count_current();
    malloc_count_reset_peak();
#endif
#endif

    // built rmq
    m_rmq = lce::rmq::rmq_n<t_index_type>(m_lcp);

#ifdef LCE_BENCHMARK_INTERNAL
    fmt::print(" rmq_time={}", t.get_and_reset());
#ifdef LCE_BENCHMARK_SPACE
    fmt::print(" rmq_mem={}", malloc_count_current() - mem_before);
    fmt::print(" rmq_mem_peak={}", malloc_count_peak() - mem_before);
#endif
#endif
  }

  template <typename C>
  lce_classic(C const& container)
      : lce_classic(container.data(), container.size()) {
  }

  // Return the number of common letters in text[i..] and text[j..].
  size_t lce(size_t i, size_t j) const {
    if (i == j) [[unlikely]] {
      assert(i < m_size);
      return m_size - i;
    }
    return lce_uneq(i, j);
  }

  // Return the number of common letters in text[i..] and text[j..]. Here i and
  // j must be different.
  size_t lce_uneq(size_t i, size_t j) const {
    assert(i != j);
    return lce_lr(i, j);
  }

  // Return the number of common letters in text[i..] and text[j..].
  // Here l must be smaller than r.
  size_t lce_lr(size_t l, size_t r) const {
    return m_lcp[m_rmq.rmq_shifted(m_isa[l], m_isa[r])];
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

 private:
  std::vector<t_index_type> m_isa;
  std::vector<t_index_type> m_lcp;

  char_type const* m_text;
  size_t m_size;
  lce::rmq::rmq_n<t_index_type> m_rmq;
};
}  // namespace lce::ds
