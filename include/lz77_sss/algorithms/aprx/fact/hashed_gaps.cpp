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

template <typename text_t>
inline lz77_sss::factor lz77_sss::factorizer<text_t>::longest_prev_occ(
    const par_buf_t& buf, uint64_t& gap, uint64_t pos)
{
    factor f = factor::literal(T.to_char(T[pos]));

    while (gap < buf.gaps.size() && buf.gaps[gap].end <= pos) {
        gap++;
    }

    if (gap == buf.gaps.size() || buf.gaps[gap].beg > pos) {
        return f;
    }

    const uint32_t* refs = buf.slot_refs.data() + buf.gaps[gap].ref_beg + (pos - buf.gaps[gap].beg);

    uint64_t prev_src = pos;

    for (int64_t x = num_patt_lens - 1; x >= 0; x--) {
        const uint32_t k = refs[x * buf.num_pos];
        if (k == rh_idx_t::no_slot) continue;
        const uint64_t src = buf.vals[k];

        if (src < pos && src != prev_src && T[src] == T[pos]) {
            prev_src = src;
            const uint64_t len = LCE_R(src, pos);

            if (len > f.len || (len == f.len && src > f.src)) {
                f.src = src;
                f.len = len;
            }
        }
    }

    for (const uint64_t dist : { buf.gaps[gap].dist_prev, buf.gaps[gap].dist_next }) {
        if (dist == 0 || dist > pos || pos - dist == prev_src || T[pos - dist] != T[pos]) continue;
        const uint64_t len = LCE_R(pos - dist, pos);

        if (len > f.len || (len == f.len && pos - dist > f.src)) {
            f.src = pos - dist;
            f.len = len;
        }
    }

    #ifndef NDEBUG
    for (uint64_t j = 0; j < f.len; j++) {
        assert(T[pos + j] == T[f.src + j]);
    }
    #endif

    return f;
}

template <typename text_t>
template <typename next_lpf_t>
void lz77_sss::factorizer<text_t>::prepare_block(
    next_lpf_t next_lpf, const par_blk_t& blk, uint64_t end,
    par_buf_t& buf, uint64_t num_owners, std::vector<uint64_t>& prefix_fps)
{
    buf.gaps.clear();
    lpf_cursor_t lpf_cursor = blk.lpf_cursor;
    lpf_phrase phr = next_lpf(lpf_cursor);
    uint64_t num_pos = 0;
    uint64_t dist_prev = blk.dist_prev;

    for (uint64_t i = blk.beg; i < end;) {
        const uint64_t gap_end = std::min<uint64_t>(phr.beg, end);
        const uint64_t dist_next = phr.beg < n ? uint64_t(phr.beg) - uint64_t(phr.src) : 0;

        if (i < gap_end) {
            const uint64_t gap_range_end = phr.beg < end ? gap_end + 1 : gap_end;
            buf.gaps.emplace_back(par_gap_t { .beg = i, .end = gap_range_end, .ref_beg = num_pos,
                .dist_prev = dist_prev, .dist_next = dist_next });
            num_pos += gap_range_end - i;
        }

        if (phr.beg >= end) break;
        i = phr.end;
        dist_prev = dist_next;
        phr = next_lpf(lpf_cursor);
    }

    buf.num_pos = num_pos;
    no_init_resize(buf.slot_refs, num_pos * num_patt_lens);

    for (const par_gap_t& gap : buf.gaps) {
        rh_idx.fill_slots(gap.beg, gap.end, buf.slot_refs.data() + gap.ref_beg, num_pos, prefix_fps);
    }

    if (num_owners == 0) {
        exchange_block(buf);
        return;
    }

    buf.owner_beg.assign(num_owners + 1, 0);

    for (uint64_t k = 0; k < num_pos * num_patt_lens; k++) {
        if (buf.slot_refs[k] != rh_idx_t::no_slot) {
            buf.owner_beg[rh_idx.owner(buf.slot_refs[k]) + 1]++;
        }
    }

    for (uint64_t o = 1; o <= num_owners; o++) {
        buf.owner_beg[o] += buf.owner_beg[o - 1];
    }

    buf.owner_fill.assign(buf.owner_beg.begin(), buf.owner_beg.end() - 1);
    no_init_resize(buf.slots, buf.owner_beg[num_owners]);
    no_init_resize(buf.vals, buf.owner_beg[num_owners]);

    for (const par_gap_t& gap : buf.gaps) {
        for (uint64_t pos = gap.beg; pos < gap.end; pos++) {
            uint32_t* refs = buf.slot_refs.data() + gap.ref_beg + (pos - gap.beg);

            for (int64_t x = num_patt_lens - 1; x >= 0; x--) {
                const uint32_t s = refs[x * num_pos];

                if (s != rh_idx_t::no_slot) {
                    const uint64_t dst = buf.owner_fill[rh_idx.owner(s)]++;
                    buf.slots[dst] = s;
                    buf.vals[dst] = pos;
                    refs[x * num_pos] = dst;
                }
            }
        }
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::exchange_block(par_buf_t& buf)
{
    const uint64_t num_pos = buf.num_pos;
    no_init_resize(buf.vals, num_pos * num_patt_lens);
    uint64_t dst = 0;

    for (const par_gap_t& gap : buf.gaps) {
        for (uint64_t pos = gap.beg; pos < gap.end; pos++) {
            const uint64_t ref = gap.ref_beg + (pos - gap.beg);
            uint32_t* refs = buf.slot_refs.data() + ref;

            if (ref + 4 < num_pos) {
                for (uint64_t x = 0; x < num_patt_lens; x++) {
                    if (refs[x * num_pos + 4] != rh_idx_t::no_slot) {
                        rh_idx.prefetch(refs[x * num_pos + 4]);
                    }
                }
            }

            for (int64_t x = num_patt_lens - 1; x >= 0; x--) {
                const uint32_t s = refs[x * num_pos];

                if (s != rh_idx_t::no_slot) {
                    buf.vals[dst] = rh_idx.exchange(s, pos);
                    refs[x * num_pos] = dst++;
                }
            }
        }
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::update_rh_idx(
    std::vector<par_buf_t>& bufs, uint64_t owner, uint64_t num_blks_in_round)
{
    for (uint64_t b = 0; b < num_blks_in_round; b++) {
        par_buf_t& buf = bufs[b];
        const uint64_t beg = buf.owner_beg[owner];
        const uint64_t end = buf.owner_beg[owner + 1];

        for (uint64_t k = beg; k < end; k++) {
            if (k + 16 < end) rh_idx.prefetch(buf.slots[k + 16]);
            buf.vals[k] = rh_idx.exchange(buf.slots[k], buf.vals[k]);
        }
    }
}

template <typename text_t>
void lz77_sss::factorizer<text_t>::output_blocks(
    merging_output& output, const std::vector<par_blk_t>& blks,
    const std::vector<std::array<std::vector<factor>, par_fact_bufs>>& facts,
    uint64_t round, uint64_t blks_per_round, uint64_t num_blks_in_round, uint64_t& out_end)
{
    for (uint64_t b = 0; b < num_blks_in_round; b++) {
        uint64_t pos = blks[round * blks_per_round + b].beg;

        for (factor f : facts[b][round % par_fact_bufs]) {
            uint64_t len = f.text_len();

            if (pos + len <= out_end) {
                pos += len;
                continue;
            }

            if (pos < out_end) {
                const uint64_t cut_left = out_end - pos;
                f.src = f.src + cut_left;
                f.len = f.len - cut_left;
                len -= cut_left;
                pos = out_end;
            }

            #ifndef NDEBUG
            assert((f.is_literal() && f.src == T.to_char(T[pos])) || f.src < pos);
            assert(f.len <= n - pos);

            for (uint64_t x = 0; x < f.len; x++) {
                assert(T[f.src + x] == T[pos + x]);
            }
            #endif

            output.add(pos, f);
            pos += len;
            out_end = pos;
        }
    }
}

template <typename text_t>
template <typename next_lpf_t>
void lz77_sss::factorizer<text_t>::factorize_block(
    next_lpf_t next_lpf, const par_blk_t& blk, uint64_t end,
    const par_buf_t& buf, std::vector<factor>& facts)
{
    facts.clear();
    lpf_cursor_t lpf_cursor = blk.lpf_cursor;
    lpf_phrase phr = next_lpf(lpf_cursor);
    uint64_t gap = 0;

    for (uint64_t i = blk.beg; true;) {
        uint64_t gap_end = phr.beg;

        if (i < gap_end) {
            do {
                factor f = longest_prev_occ(buf, gap, i);
                i += f.text_len();

                if (i > gap_end) {
                    if (i <= phr.end) {
                        f.len = f.len - (i - gap_end);
                        i = gap_end;
                    } else {
                        do {
                            phr = next_lpf(lpf_cursor);
                        } while (phr.end <= i);

                        gap_end = phr.beg;
                    }
                }

                facts.emplace_back(f);
            } while (i < gap_end && i < end);
        }

        if (i >= end) {
            break;
        }

        uint64_t cut_left = i - gap_end;

        factor best {
            .src = phr.src + cut_left,
            .len = (phr.end - phr.beg) - cut_left
        };

        factor f = longest_prev_occ(buf, gap, i);

        if (f.len > best.len) {
            best = f;
        }

        #ifndef NDEBUG
        assert(best.src < i);

        for (uint64_t j = 0; j < best.len; j++) {
            assert(T[i + j] == T[best.src + j]);
        }
        #endif

        facts.emplace_back(best);
        i += best.len;

        while (phr.end <= i) {
            phr = next_lpf(lpf_cursor);
        }
    }
}

template <typename text_t>
template <typename lpf_beg_t, typename next_lpf_t>
void lz77_sss::factorizer<text_t>::factorize_hashed_gaps(
    factor_sink& output,
    lpf_beg_t lpf_beg,
    next_lpf_t next_lpf)
{
    log_phase_begin(log, "factorizing gaps (approximate)");

    std::vector<par_blk_t> blks;
    uint64_t max_blk_num_pos = 0;

    {
        uint64_t max_blk_gaps = 0;
        lpf_cursor_t lpf_cursor = lpf_beg();
        uint64_t gap_beg = 0;
        uint64_t gap_offs = 0;
        uint64_t nxt_blk_gap_offs = 0;
        uint64_t blk_gaps = 0;
        uint64_t dist_prev = 0;

        while (true) {
            const lpf_cursor_t phr_cursor = lpf_cursor;
            const lpf_phrase phr = next_lpf(lpf_cursor);
            const uint64_t gap_len = phr.beg > gap_beg ? phr.beg - gap_beg : 0;

            while (nxt_blk_gap_offs < gap_offs + gap_len) {
                max_blk_gaps = std::max(max_blk_gaps, blk_gaps);
                blk_gaps = 0;

                blks.emplace_back(par_blk_t {
                    .beg = gap_beg + (nxt_blk_gap_offs - gap_offs),
                    .lpf_cursor = phr_cursor,
                    .dist_prev = dist_prev });

                nxt_blk_gap_offs += par_gap_blk_len;
            }

            gap_offs += gap_len;
            if (phr.beg >= n) break;
            if (gap_len > 0) blk_gaps++;
            gap_beg = std::max<uint64_t>(gap_beg, phr.end);
            dist_prev = uint64_t(phr.beg) - uint64_t(phr.src);
        }

        max_blk_num_pos = par_gap_blk_len + std::max(max_blk_gaps, blk_gaps);
    }

    const uint64_t num_blks = blks.size();
    blks.emplace_back(par_blk_t { .beg = n, .lpf_cursor = lpf_cursor_t { }, .dist_prev = 0 });

    const bool has_out_thr = p >= 3;
    const uint64_t num_workers = has_out_thr ? p - 1 : p;
    const uint64_t blks_per_round = par_blks_per_thr * num_workers;
    const uint64_t num_owners = num_workers == 1 ? 0 : par_owners_per_thr * num_workers;
    const uint64_t num_rounds = div_ceil<uint64_t>(num_blks, blks_per_round);
    if (num_owners > 0) rh_idx.set_num_owners(num_owners);

    auto blks_in_round = [&](uint64_t round) {
        return std::min<uint64_t>(blks_per_round, num_blks - round * blks_per_round);
    };

    std::vector<par_buf_t> bufs(blks_per_round);
    std::vector<std::array<std::vector<factor>, par_fact_bufs>> facts(blks_per_round);
    spin_barrier barrier(num_workers);
    std::vector<std::atomic<uint8_t>> claimed(std::max<uint64_t>(blks_per_round, num_owners));
    std::atomic<uint64_t> rounds_factorized = 0;
    std::atomic<uint64_t> rounds_output = 0;
    uint64_t out_end = 0;
    merging_output merged(output, num_fact);

    auto claim = [&](uint64_t t) {
        return claimed[t].load(std::memory_order_relaxed) == 0 &&
            claimed[t].exchange(1, std::memory_order_relaxed) == 0;
    };

    auto run_tasks = [&](uint64_t worker, uint64_t num_tasks, auto fnc) {
        for (uint64_t t = worker; t < num_tasks; t += num_workers) {
            if (claim(t)) fnc(t);
        }

        for (uint64_t i = 1; i < num_tasks; i++) {
            const uint64_t t = (worker + i) % num_tasks;
            if (claim(t)) fnc(t);
        }
    };

    auto reset_tasks = [&]() {
        for (std::atomic<uint8_t>& c : claimed) c.store(0, std::memory_order_relaxed);
    };

    phase_progress progress(log, num_rounds);

    #pragma omp parallel num_threads(p)
    {
        const uint16_t i_p = omp_get_thread_num();

        if (has_out_thr && i_p == 0) {
            for (uint64_t round = 0; round < num_rounds; round++) {
                spin_until([&]() { return rounds_factorized.load(std::memory_order_acquire) > round; });
                output_blocks(merged, blks, facts, round, blks_per_round, blks_in_round(round), out_end);
                rounds_output.store(round + 1, std::memory_order_release);
                progress.reached(round + 1);
            }
        } else {
            const uint64_t worker = has_out_thr ? i_p - 1 : i_p;
            std::vector<uint64_t> prefix_fps;

            for (uint64_t b = worker; b < std::min(blks_per_round, num_blks); b += num_workers) {
                bufs[b].slot_refs.reserve(max_blk_num_pos * num_patt_lens);
                bufs[b].slots.reserve(max_blk_num_pos * num_patt_lens);
                bufs[b].vals.reserve(max_blk_num_pos * num_patt_lens);
            }

            barrier.wait();

            run_tasks(worker, blks_in_round(0), [&](uint64_t b) {
                prepare_block(next_lpf, blks[b], blks[b + 1].beg, bufs[b], num_owners, prefix_fps);
            });

            barrier.wait(reset_tasks);

            run_tasks(worker, num_owners, [&](uint64_t o) {
                update_rh_idx(bufs, o, blks_in_round(0));
            });

            barrier.wait(reset_tasks);

            for (uint64_t round = 0; round < num_rounds; round++) {
                if (has_out_thr && round >= par_fact_bufs) {
                    spin_until([&]() { return rounds_output.load(std::memory_order_acquire) + par_fact_bufs > round; });
                }

                run_tasks(worker, blks_in_round(round), [&](uint64_t b) {
                    const uint64_t blk = round * blks_per_round + b;
                    factorize_block(next_lpf, blks[blk], blks[blk + 1].beg,
                        bufs[b], facts[b][round % par_fact_bufs]);

                    if (round + 1 < num_rounds && b < blks_in_round(round + 1)) {
                        prepare_block(next_lpf, blks[blk + blks_per_round], blks[blk + blks_per_round + 1].beg,
                            bufs[b], num_owners, prefix_fps);
                    }
                });

                barrier.wait([&]() {
                    reset_tasks();
                    if (has_out_thr) rounds_factorized.store(round + 1, std::memory_order_release);
                });

                if (!has_out_thr && i_p == 0) {
                    output_blocks(merged, blks, facts, round, blks_per_round, blks_in_round(round), out_end);
                    progress.reached(round + 1);
                }

                if (round + 1 < num_rounds) {
                    run_tasks(worker, num_owners, [&](uint64_t o) {
                        update_rh_idx(bufs, o, blks_in_round(round + 1));
                    });
                }

                barrier.wait(reset_tasks);
            }
        }
    }

    merged.flush();
    log_progress_done(log);

    if (log) {
        record_phase_time("factorize", time_diff_ns(time, now()));
        time = log_runtime(time);
    }
}
