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

#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau, typename char_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::build_PSV_NSV_S()
{
    if (log) {
        std::cout << "building NSV_S and PSV_S" << std::flush;
    }

    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SA_S = LCE.get_sa_s();
    pos_t s = S.size();

    #ifndef NDEBUG
    for (uint32_t i = 1; i < s; i++) {
        assert(S[i - 1] < S[i]);
    }
    #endif

    no_init_resize(PSV_S, s);
    no_init_resize(NSV_S, s);

    if (s == 0) [[unlikely]] {
        if (log) {
            log_phase("nsv_psv", time_diff_ns(time, now()));
            time = log_runtime(time);
        }
        return;
    }

    PSV_S[0] = s;
    NSV_S[s - 1] = s;

    for (uint32_t i = 1; i < s; i++) {
        uint32_t j = i - 1;

        while (j != s && SA_S[j] > SA_S[i]) {
            NSV_S[j] = i;
            j = PSV_S[j];
        }

        PSV_S[i] = j;

        if (j != s) {
            NSV_S[j] = s;
        }

        #ifndef NDEBUG
        assert(PSV_S[i] == s || S[SA_S[PSV_S[i]]] < S[SA_S[i]]);
        #endif
    }

    #ifndef NDEBUG
    for (int64_t i = s - 1; i >= 0; i--) {
        assert(NSV_S[i] == s || S[SA_S[NSV_S[i]]] < S[SA_S[i]]);
    }
    #endif

    if (log) {
        log_phase("nsv_psv", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau, typename char_t>
void lz77_sss<pos_t>::factorizer<tau, char_t>::build_PGV_NGV_S()
{
    if (log) {
        std::cout << "building PGV_S and NGV_S" << std::flush;
    }

    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SA_S = LCE.get_sa_s();
    [[maybe_unused]] const std::vector<uint32_t>& ISA_S = LCE.get_isa_s();
    pos_t s = S.size();

    #ifndef NDEBUG
    for (uint32_t i = 1; i < s; i++) {
        assert(S[i - 1] < S[i]);
    }
    #endif

    no_init_resize(PGV_S, s);
    no_init_resize(NGV_S, s);

    if (s == 0) [[unlikely]] {
        if (log) {
            log_phase("ngv_pgv", time_diff_ns(time, now()));
            time = log_runtime(time);
        }
        return;
    }

    PGV_S[0] = s;
    NGV_S[s - 1] = s;

    for (uint32_t i = 1; i < s; i++) {
        uint32_t j = i - 1;

        while (j != s && SA_S[j] < SA_S[i]) {
            NGV_S[j] = i;
            j = PGV_S[j];
        }

        PGV_S[i] = j;
        
        if (j != s) {
            NGV_S[j] = s;
        }

        #ifndef NDEBUG
        assert(PGV_S[i] == s || S[SA_S[PGV_S[i]]] > S[SA_S[i]]);
        #endif
    }

    #ifndef NDEBUG
    for (int64_t i = s - 1; i >= 0; i--) {
        assert(NGV_S[i] == s || S[SA_S[NGV_S[i]]] > S[SA_S[i]]);
    }
    #endif

    if (log) {
        log_phase("ngv_pgv", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}