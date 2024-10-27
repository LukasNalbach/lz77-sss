#pragma once

#include <lz77_sss/lz77_sss.hpp>

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::build_PSV_NSV_S()
{
    if (log) {
        std::cout << "building NSV_S and PSV_S" << std::flush;
    }

    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SSA_S = LCE.get_ssa();
    pos_t s = S.size();

    #ifndef NDEBUG
    for (uint32_t i = 1; i < s; i++) {
        assert(S[i - 1] < S[i]);
    }
    #endif

    no_init_resize(PSV_S, s);
    no_init_resize(NSV_S, s);
    PSV_S[0] = 0;
    NSV_S[s - 1] = s;

    for (uint32_t i = 1; i < s; i++) {
        uint32_t j = i - 1;

        while (j != 0 && SSA_S[j] > SSA_S[i]) {
            NSV_S[j] = i;
            j = PSV_S[j];
        }

        PSV_S[i] = j;
        NSV_S[j] = s;

        #ifndef NDEBUG
        assert(PSV_S[i] == 0 || S[SSA_S[PSV_S[i]]] < S[SSA_S[i]]);
        #endif
    }

    #ifndef NDEBUG
    for (int64_t i = s - 1; i >= 0; i--) {
        assert(NSV_S[i] == s || S[SSA_S[NSV_S[i]]] < S[SSA_S[i]]);
    }
    #endif

    if (log) {
        log_phase("nsv_psv", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::build_PGV_NGV_S()
{
    if (log) {
        std::cout << "building PGV_S and NGV_S" << std::flush;
    }

    const std::vector<pos_t>& S = LCE.get_sync_set();
    const std::vector<uint32_t>& SSA_S = LCE.get_ssa();
    const std::vector<uint32_t>& ISSA_S = LCE.get_issa();
    pos_t s = S.size();

    #ifndef NDEBUG
    for (uint32_t i = 1; i < s; i++) {
        assert(S[i - 1] < S[i]);
    }
    #endif

    no_init_resize(PGV_S, s);
    no_init_resize(NGV_S, s);
    PGV_S[0] = 0;
    NGV_S[s - 1] = s;

    for (uint32_t i = 1; i < s; i++) {
        uint32_t j = i - 1;

        while (j != 0 && SSA_S[j] < SSA_S[i]) {
            NGV_S[j] = i;
            j = PGV_S[j];
        }

        PGV_S[i] = j;
        NGV_S[j] = s;

        #ifndef NDEBUG
        assert(PGV_S[i] == 0 || S[SSA_S[PGV_S[i]]] > S[SSA_S[i]]);
        #endif
    }

    #ifndef NDEBUG
    for (int64_t i = s - 1; i >= 0; i--) {
        assert(NGV_S[i] == s || S[SSA_S[NGV_S[i]]] > S[SSA_S[i]]);
    }
    #endif

    if (log) {
        log_phase("ngv_pgv", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}