/**
 * lz77/gzip9_factorizer.hpp
 * part of pdinklag/lz77
 * 
 * MIT License
 * 
 * Copyright (c) Patrick Dinklage
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

#ifndef _LZ77_GZIP9_FACTOR_HPP
#define _LZ77_GZIP9_FACTOR_HPP

#include <algorithm>
#include <cassert>
#include <concepts>
#include <iterator>
#include <cstring>

#include "factor.hpp"

namespace lz77 {

/**
 * \brief Computes the `gzip -9` factorization of the input
 * 
 * This implementation produces exactly the same factorization as `gzip -9` without any subsequent encoding steps.
 * It uses a sliding window and a number of heuristics to ensure a small encoding in the deflate format.
 */
class Gzip9Factorizer {
private:
    static constexpr bool track_stats_ = false;

    struct Stats {
        size_t chain_num;
        size_t chain_length_max;
        size_t chain_length_sum;
        size_t nice_num;
        size_t good_num;
        size_t greedy_skips;
        size_t match_ops;
    };

    // the following parameters pertain to gzip -9
    // other configurations may be implemented in the future
    static constexpr size_t min_match_ = 3;
    static constexpr size_t max_match_ = 258;
    static constexpr size_t window_bits_ = 15;
    static constexpr size_t window_size_ = 1ULL << window_bits_;
    static constexpr size_t buf_capacity_ = 2 * window_size_;
    static constexpr size_t window_mask_ = window_size_ - 1;
    static constexpr size_t num_chains_ = 1ULL << window_bits_;
    static constexpr size_t chain_mask_ = num_chains_ - 1;
    static constexpr size_t hash_shift_ = 5; // if configured to window_bits_ / min_match_, hash function becomes rolling (but isn't used as such)
    static constexpr size_t min_lookahead_ = max_match_ + min_match_ + 1;
    static constexpr size_t max_dist_ = window_size_ - min_lookahead_;
    static constexpr size_t max_chain_length_ = 4096; // gzip - 9
    static constexpr size_t nice_match_ = 258; // gzip - 9
    static constexpr size_t lazy_match_ = 258; // gzip - 9
    static constexpr size_t good_match_ = 32; // gzip - 9
    static constexpr size_t good_laziness_ = 4;
    static constexpr size_t too_far_ = 4096;

    using WindowIndex = std::conditional<window_bits_ <= 15, uint16_t, uint32_t>::type;
    static constexpr auto NIL = 0ULL;

    inline size_t hash(const size_t p) const {
        size_t h = 0;
        for(size_t i = 0; i < min_match_; i++) {
            h = ((h << hash_shift_) ^ buf_[p + i]) & chain_mask_;
        }
        return h;
    }

    uint8_t* buf_;
    uintmax_t buf_offs_;  // text position of first entry in buffer
    size_t buf_avail_; // available bytes in buffer
    size_t buf_pos_; // current read position in buffer

    uintmax_t pos_;         // next text position to encode
    size_t hash_only_;   // number of positions to skip - but still hash - after emitting a reference

    size_t prev_length_;
    uintmax_t prev_src_;
    bool prev_match_exists_;

    size_t match_length_;
    uintmax_t match_src_;

    WindowIndex* hashtable_; // memory
    WindowIndex* head_; // head of hash chains
    WindowIndex* prev_; // chains

    // stats
    Stats stats_;

    template<typename CharInput>
    size_t advance(CharInput& begin, CharInput const& end, uint8_t* buf, size_t max) {
        size_t num = 0;
        while(num < max && begin != end) {
            *buf++ = (uint8_t)*begin++;
            ++num;
        }
        return num;
    }

    template<typename FactorOutput>
    void process(FactorOutput& out) {
        const size_t relative_pos = pos_ - buf_offs_;

        // insert current string
        uintmax_t src;
        {
            const auto h = hash(buf_pos_);
            src = head_[h];
            prev_[relative_pos & window_mask_] = src;

            assert(relative_pos < buf_capacity_);
            head_[h] = relative_pos;
        }
        
        if(hash_only_) {
            --hash_only_;
        } else {
            // store previous match
            prev_length_ = match_length_;
            prev_src_ = match_src_;
            match_length_ = min_match_ - 1; // init to horrible

            // find the longest match
            if(src != NIL && prev_length_ < lazy_match_ && relative_pos - src <= max_dist_) {
                {
                    const auto limit = relative_pos > max_dist_ ? relative_pos - max_dist_ : NIL;

                    if constexpr(track_stats_) {
                        if(prev_length_ >= good_match_) ++stats_.good_num;
                    }

                    size_t chain = (prev_length_ >= good_match_) ? (max_chain_length_ / good_laziness_) : max_chain_length_;
                    size_t chain_length = 0; // for stats only

                    match_length_ = prev_length_; // we want to beat the previous match at least

                    const uint8_t* const match_begin = buf_ + buf_pos_;
                    const uint8_t* const match_end = ((buf_pos_ + max_match_ <= buf_avail_) ? match_begin + max_match_ : buf_ + buf_avail_);
                    const uint16_t prefix = *(const uint16_t*)match_begin;
                    uint16_t suffix = *(const uint16_t*)(match_begin + match_length_ - 1);

                    do {
                        if constexpr(track_stats_) ++stats_.match_ops;

                        // prepare match
                        const uint8_t* p = match_begin;
                        const uint8_t* q = buf_ + src;
                        assert(q < p);

                        // if first two characters don't match OR we cannot become better, then don't even bother
                        if(prefix == *(const uint16_t*)q && suffix == *(const uint16_t*)(q + match_length_ - 1)) {
                            // already matched first two, so skipping the first two by using the += operator below is safe
                            // the next bytes up to min_match_ must also match because we are in the corresponding hash chain
                            ++p;
                            ++q;

                            // there are at most 256 bytes left to match
                            // we only test for crossing the boundary every 8 bytes
                            static_assert((max_match_ - 2) % 8 == 0, "the funny tricks only work for max_match_ = 8k + 2 for some k");

                            while(
                                *(const uint16_t*)(p+=2) == *(const uint16_t*)(q+=2) &&
                                *(const uint16_t*)(p+=2) == *(const uint16_t*)(q+=2) &&
                                *(const uint16_t*)(p+=2) == *(const uint16_t*)(q+=2) &&
                                *(const uint16_t*)(p+=2) == *(const uint16_t*)(q+=2) &&
                                p + 1 < match_end
                            );

                            if(p < match_end && *p == *q) ++p; // final comparison

                            //const size_t length = std::min((size_t)(p - match_begin), max_match_);
                            const size_t length = (size_t)(p - match_begin);

                            // check match
                            if(length > match_length_) {
                                match_src_ = src;
                                match_length_ = length;

                                if(length >= nice_match_) {
                                    // immediately break when finding a nice match
                                    if constexpr(track_stats_) ++stats_.nice_num;
                                    break;
                                }

                                suffix = *(const uint16_t*)(buf_ + buf_pos_ + match_length_ - 1);
                            }
                        }

                        // advance in chain
                        if constexpr(track_stats_) ++chain_length;
                    } while(--chain && (src = prev_[src & window_mask_]) > limit);

                    // stats
                    if constexpr(track_stats_) {
                        stats_.chain_length_max = std::max(stats_.chain_length_max, chain_length);
                        stats_.chain_length_sum += chain_length;
                        ++stats_.chain_num;
                    }
                }

                // make match source global
                match_src_ += buf_offs_;

                // ignore a minimum match if it is too distant
                if(match_length_ == min_match_ && pos_ - match_src_ > too_far_) {
                    --match_length_;
                }
            }

            // compare current match against previous match
            if(prev_length_ >= min_match_ && match_length_ <= prev_length_) {
                // previous match was better than current, emit
                *out++ = Factor(pos_ - 1 - prev_src_, prev_length_);
                hash_only_ = prev_length_ - 2; // current and previous positions are already hashed

                // reset
                match_length_ = min_match_ - 1;
                prev_match_exists_ = false;
            } else if(prev_match_exists_) {
                // current match is better, truncate previous match to a single literal
                assert(buf_pos_ > 0);

                if constexpr(track_stats_) {
                    if(prev_length_ >= min_match_) ++stats_.greedy_skips;
                }

                *out++ = Factor((char)buf_[buf_pos_ - 1]);
            } else {
                // there is no previous match to compare with, wait for next step
                prev_match_exists_ = true;
            }
        }
    }

public:
    Gzip9Factorizer() {
        const size_t bufsize = buf_capacity_ + min_lookahead_;
        buf_ = new uint8_t[bufsize];
        for(size_t i = 0; i < bufsize; i++) buf_[i] = 0;
        
        hashtable_ = new WindowIndex[num_chains_ + window_size_];

        head_ = hashtable_;
        for(size_t i = 0; i < num_chains_; i++) head_[i] = NIL;
        
        prev_ = hashtable_ + num_chains_;
        for(size_t i = 0; i < window_size_; i++) prev_[i] = NIL;
    }

    ~Gzip9Factorizer() {
        delete[] buf_;
        delete[] hashtable_;
    }

    template<std::input_iterator Input, std::output_iterator<Factor> Output>
    requires (sizeof(std::iter_value_t<Input>) == 1)
    void factorize(Input begin, Input const& end, Output out) {
        if constexpr(track_stats_) {
            stats_.chain_length_max = 0;
            stats_.chain_length_sum = 0;
            stats_.chain_num = 0;
            stats_.nice_num = 0;
            stats_.good_num = 0;
            stats_.greedy_skips = 0;
            stats_.match_ops = 0;
        }

        // initialize
        buf_offs_ = 0;
        buf_avail_ = 0;
        buf_pos_ = 0;

        match_src_ = NIL;
        match_length_ = 0;

        prev_src_ = NIL;
        prev_length_ = 0;
        prev_match_exists_ = false;

        hash_only_ = 0;
        
        pos_ = 0;

        // fill buffer
        buf_avail_ = advance(begin, end, buf_, buf_capacity_);
        while(begin != end) {
            // process while buffer has enough bytes left
            assert(buf_avail_ > min_lookahead_);
            {
                const size_t buf_border = buf_avail_ - min_lookahead_;
                while(buf_pos_ < buf_border) {
                    process(out);
                    
                    ++buf_pos_;
                    ++pos_;
                }
            }

            // buffer ran short of min lookahead, slide
            {
                assert(buf_pos_ >= window_size_);

                std::memcpy(buf_, buf_ + window_size_, window_size_);
                buf_pos_ -= window_size_;
                buf_offs_ += window_size_;

                // read more
                const size_t num_read = advance(begin, end, buf_ + window_size_, window_size_);
                buf_avail_ = window_size_ + num_read;
            }

            // clean up hash chains
            {
                for(size_t i = 0; i < num_chains_; i++) {
                    const auto m = head_[i];
                    head_[i] = (m >= window_size_) ? m - window_size_ : NIL;
                }

                for(size_t i = 0; i < window_size_; i++) {
                    const auto m = prev_[i];
                    prev_[i] = (m >= window_size_) ? m - window_size_ : NIL;
                }
            }
        }

        // process final window
        while(buf_pos_ + min_match_ <= buf_avail_) {
            process(out);
            
            ++pos_;
            ++buf_pos_;
        }

        // emit remaining literals
        while(buf_pos_ < buf_avail_) {
            if(hash_only_) {
                --hash_only_;
            } else {
                *out++ = Factor((char)buf_[buf_pos_]);
            }

            ++buf_pos_;
            ++pos_;
        }
    }
};

}

#endif
