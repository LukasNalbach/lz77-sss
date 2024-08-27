#pragma once

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, out_it_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::build_c(std::istream_iterator<factor>& ifile_approx_it) {
    if (log) {
        std::cout << "setting delta = " << delta << std::endl;
        std::cout << "building C" << std::flush;
    }

    C.reserve(num_phr + n / delta);
    C.emplace_back(0);
    pos_t end_last_phr = 0;
    pos_t end_next_phr;

    for (pos_t i = 1; i < num_phr; i++) {
        end_next_phr = end_last_phr + std::max<pos_t>(1, ifile_approx_it++->len);

        while (end_next_phr - end_last_phr > delta) {
            end_last_phr += delta;
            C.emplace_back(end_last_phr);
        }

        C.emplace_back(end_next_phr);
        end_last_phr = end_next_phr;
    }

    C.shrink_to_fit();
    c = C.size();

    if (log) {
        std::cout << " (" << format_size(c * sizeof(pos_t)) << ")";
        time = log_runtime(time);
        std::cout << "c / num. of phrases = " << c / (double) num_phr << std::endl;
    }
}

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, out_it_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::build_idx_C() {
    if (log) {
        std::cout << "building sample-index for C:" << std::endl;
    }

    idx_C.build(T, C, LCE, delta, n / (1.0 * num_phr), transf_mode == optimized, log);

    if (log) {
        std::cout << "size: " << format_size(idx_C.size_in_bytes());
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, out_it_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::build_p() {
    if (log) {
        std::cout << "computing P" << std::flush;
    }

    P.reserve(c);

    for (sidx_t i = 0; i < c; i++) {
        if constexpr (is_static<range_ds_t>()) {
            P.emplace_back(point_t{.weight = i});
        } else {
            P.emplace_back(point_t{});
        }
    }

    for (sidx_t i = 0; i < c; i++) {
        P[idx_C.spa(i)].x = i;
        P[idx_C.ssa(i)].y = i;
    }

    if (log) {
        std::cout << " (" << format_size(sizeof(point_t) * c) << ")";
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, out_it_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::build_ps_sp() {
    if (log) {
        std::cout << "computing PS and SP" << std::flush;
    }
    
    PS.reserve(c);

    for (sidx_t i = 0; i < c; i++) {
        PS.push_back(P[idx_C.spa(i)].y);
    }

    SP.reserve(c);

    for (sidx_t i = 0; i < c; i++) {
        SP.push_back(P[idx_C.ssa(i)].x);
    }

    if (log) {
        std::cout << " (" << format_size(2 * c * sizeof(sidx_t)) << ")";
        time = log_runtime(time);
    }
}

template <typename pos_t>
template <uint64_t tau, typename out_it_t>
template <typename sidx_t, transform_mode transf_mode, template <typename> typename range_ds_t>
void lz77_sss<pos_t>::factorizer<tau, out_it_t>::exact_factorizer<sidx_t, transf_mode, range_ds_t>::adjust_sample_index(sidx_t& idx, pos_t pos) {
    while (idx < c && C[idx] < pos) {
        idx++;
    }

    while (idx > 0 && C[idx - 1] >= pos) {
        idx--;
    }
}