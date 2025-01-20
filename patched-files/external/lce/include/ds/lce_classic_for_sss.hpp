/*******************************************************************************
 * lce/ds/lce_classic_for_sss.hpp
 *
 * Copyright (C) 2022 Alexander Herlez <alexander.herlez@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#include <assert.h>

#include <cstdint>
#include <gsaca-double-sort-par.hpp>

#include "ds/lce_naive_wordwise.hpp"
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

template <typename t_index_type, size_t t_tau>
class lce_classic_for_sss {
 public:
  lce_classic_for_sss() : m_size{0} {
  }

  lce_classic_for_sss(uint8_t const* text, size_t text_size,
                      t_index_type const* reduced_fps, size_t reduced_fps_size,
                      std::vector<t_index_type> const& sss)
      : m_size(reduced_fps_size) {
    m_sa.resize(reduced_fps_size);
    // sort sa
#ifdef LCE_BENCHMARK_INTERNAL
    lce::util::timer t;
#ifdef LCE_BENCHMARK_SPACE
    size_t mem_before = malloc_count_current();
    malloc_count_reset_peak();
#endif
#endif
    gsaca_for_lce(reduced_fps, m_sa.data(), reduced_fps_size);

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
    m_isa.resize(reduced_fps_size);
#pragma omp parallel for
    for (size_t i = 0; i < m_sa.size(); ++i) {
      m_isa[m_sa[i]] = i;
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
    m_lcp.resize(m_sa.size());
    m_lcp[0] = 0;
    size_t current_lcp = 0;

#pragma omp parallel
    {
      const int t = omp_get_thread_num();
      const int nt = omp_get_num_threads();
      const size_t slice_size = reduced_fps_size / nt;

      const size_t begin = t * slice_size;
      const size_t end = (t < nt - 1) ? (t + 1) * slice_size : reduced_fps_size;

      size_t current_lcp = 0;
      for (size_t i{begin}; i < end; ++i) {
        size_t suffix_array_pos = m_isa[i];
        if (suffix_array_pos == 0) {
          continue;
        }
        assert(suffix_array_pos != 0);

        size_t preceding_suffix_pos = m_sa[suffix_array_pos - 1];
        current_lcp += lce_naive_wordwise_xor<uint8_t>::lce_uneq(
            text, text_size, sss[i] + current_lcp,
            sss[preceding_suffix_pos] + current_lcp);
        m_lcp[suffix_array_pos] = current_lcp;
        assert(lce_naive_wordwise_xor<uint8_t>::lce_uneq(
                   text, text_size, sss[i], sss[preceding_suffix_pos]) ==
               current_lcp);
        if (i == end-1) break;
        uint64_t diff = sss[i + 1] - sss[i];
        if (current_lcp < 2 * t_tau + diff) {
          current_lcp = 0;
        } else {
          current_lcp -= diff;
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

    // build rmq
    m_rmq = lce::rmq::rmq_n<t_index_type>(m_lcp);

#ifdef LCE_BENCHMARK_INTERNAL
    fmt::print(" rmq_time={}", t.get_and_reset());
#ifdef LCE_BENCHMARK_SPACE
    fmt::print(" rmq_mem={}", malloc_count_current() - mem_before);
    fmt::print(" rmq_mem_peak={}", malloc_count_peak() - mem_before);
#endif
#endif
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

  const std::vector<uint32_t>& get_ssa() {
    return m_sa;
  }

  const std::vector<uint32_t>& get_issa() {
    return m_isa;
  }

  void delete_ssa() {
    m_sa.clear();
    m_sa.shrink_to_fit();
  }

 private:
  size_t m_size;
  std::vector<uint32_t> m_sa;
  std::vector<uint32_t> m_isa;
  std::vector<t_index_type> m_lcp;
  lce::rmq::rmq_n<t_index_type> m_rmq;
};
}  // namespace lce::ds
