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

#define LZ77_SSS_DECLARATIONS_ONLY
#include <lz77_sss/lz77_sss.hpp>
#include <lz77_sss/algorithms/aprx_to_exact/common.cpp>
#include <lz77_sss/algorithms/aprx_to_exact/with_interval_samples.cpp>
#include <lz77_sss/algorithms/aprx_to_exact/without_interval_samples.cpp>

template void lz77_sss::factorizer<lz77_sss::LZ77_SSS_INST_TEXT>::exact_transformer::build_C();
template void lz77_sss::factorizer<lz77_sss::LZ77_SSS_INST_TEXT>::exact_transformer::build_idx_C();
template void lz77_sss::factorizer<lz77_sss::LZ77_SSS_INST_TEXT>::exact_transformer::build_P();
template void lz77_sss::factorizer<lz77_sss::LZ77_SSS_INST_TEXT>::exact_transformer::build_Pi_Psi();
template void lz77_sss::factorizer<lz77_sss::LZ77_SSS_INST_TEXT>::exact_transformer::transform_to_exact_with_interval_samples(lz77_sss::factor_sink&);
template void lz77_sss::factorizer<lz77_sss::LZ77_SSS_INST_TEXT>::exact_transformer::transform_to_exact_without_interval_samples(lz77_sss::factor_sink&);
