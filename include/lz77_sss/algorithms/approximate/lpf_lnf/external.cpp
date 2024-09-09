#pragma once

#include <lz77_sss/lz77_sss.hpp>
#include <gtl/btree.hpp>

template <typename pos_t>
template <uint64_t tau>
void lz77_sss<pos_t>::factorizer<tau>::build_LPF_all_external() {
    struct lpf_comp {
        bool operator()(const lpf &p1, const lpf &p2) const {
            return p1.beg < p2.beg ||
              (p1.beg == p2.beg && p1.end > p2.end);
        }
    };

    gtl::btree_set<lpf, lpf_comp> T_LPF;
    std::vector<lpf> LPF_buf;
    lpf_file_name = "lpf_" + random_alphanumeric_string(10);
    std::ofstream lpf_ofile(lpf_file_name);
    pos_t tau_ = tau;

    build_LPF_all([&](lpf&& p, pos_t i){
        T_LPF.emplace_hint(T_LPF.end(), p);

        if (T_LPF.size() > max_buffered_lpf_phrases) {
            auto it_end = T_LPF.upper_bound(lpf {
                .beg = i - tau_,
                .end = n,
                .src = 0
            });

            LPF_buf.insert(LPF_buf.end(),
                T_LPF.begin(), it_end);

            T_LPF.erase(T_LPF.begin(), it_end);

            lpf_ofile.write((char*) &LPF_buf[0],
                LPF_buf.size() * sizeof(lpf));

            LPF_buf.clear();
        }
    });

    LPF_buf.insert(LPF_buf.end(),
        T_LPF.begin(), T_LPF.end());

    lpf_ofile.write((char*) &LPF_buf[0],
        LPF_buf.size() * sizeof(lpf));

    LPF_buf.clear();
    lpf_ofile.close();
}