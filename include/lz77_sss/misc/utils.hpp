#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <climits>
#include <string>
#include <functional>
#include <unistd.h>

#include <omp.h>
#include <malloc_count/malloc_count.h>
#include <rolling_hash/rolling_hash.hpp>

using uint128_t = alx::rolling_hash::uint128_t;

std::chrono::steady_clock::time_point now() {
    return std::chrono::steady_clock::now();
}

std::string format_time(uint64_t ns) {
    std::string time_str;

    if (ns > 10000000000) {
        time_str = std::to_string(ns/1000000000) + " s";
    } else if (ns > 10000000) {
        time_str = std::to_string(ns/1000000) + " ms";
    } else if (ns > 10000) {
        time_str = std::to_string(ns/1000) + " us";
    } else {
        time_str = std::to_string(ns) + " ns";
    }

    return time_str;
}

std::string format_size(uint64_t B) {
    std::string size_str;

    if (B > 10000000000) {
        size_str = std::to_string(B/1000000000) + " GB";
    } else if (B > 10000000) {
        size_str = std::to_string(B/1000000) + " MB";
    } else if (B > 10000) {
        size_str = std::to_string(B/1000) + " KB";
    } else {
        size_str = std::to_string(B) + " B";
    }
    
    return size_str;
}

double throughput(uint64_t bytes, uint64_t ns) {
    return 1'000.0 * (bytes / (double) ns);
}

std::string format_throughput(uint64_t bytes, uint64_t ns) {
    return std::to_string(throughput(bytes, ns)) + " MB/s";
}

uint64_t time_diff_min(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    return std::chrono::duration_cast<std::chrono::minutes>(t2-t1).count();
}

uint64_t time_diff_min(std::chrono::steady_clock::time_point t) {
    return time_diff_min(t,std::chrono::steady_clock::now());
}

uint64_t time_diff_ns(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
}

uint64_t time_diff_ns(std::chrono::steady_clock::time_point t) {
    return time_diff_ns(t,std::chrono::steady_clock::now());
}

std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2) {
    std::cout << ", in ~ " << format_time(time_diff_ns(t1,t2)) << std::endl;
    return std::chrono::steady_clock::now();
}

std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t) {
    return log_runtime(t,std::chrono::steady_clock::now());
}

void log_message(std::string message) {
    std::cout << message << std::flush;
}

void copy_buffered(
    std::fstream& input_stream, std::fstream& output_stream,
    std::string& buffer, uint64_t from, uint64_t to, uint64_t length, uint64_t buffer_size
) {
    if (&input_stream == &output_stream) buffer_size = std::min(buffer_size, to - from);
    buffer.resize(buffer_size);

    for (uint64_t offset = 0; offset < length;) {
        uint64_t bytes_to_copy = std::min<uint64_t>(length - offset, buffer_size);
        input_stream.seekg(from + offset, std::ios::beg);
        input_stream.read(buffer.data(), bytes_to_copy);
        output_stream.seekp(to + offset, std::ios::beg);
        output_stream.write(buffer.data(), bytes_to_copy);
        offset += bytes_to_copy;
    }
}

void read_from_file(std::istream& in, const char* data, uint64_t size) {
    uint64_t size_left = size;
    uint64_t bytes_to_read;

    while (size_left > 0) {
        bytes_to_read = std::min(size_left,(uint64_t)INT_MAX);
        in.read((char*)&data[size-size_left],bytes_to_read);
        size_left -= bytes_to_read;
    }
}

void write_to_file(std::ostream& out, const char* data, uint64_t size) {
    uint64_t size_left = size;
    uint64_t bytes_to_write;

    while (size_left > 0) {
        bytes_to_write = std::min(size_left,(uint64_t)INT_MAX);
        out.write((char*)&data[size-size_left],bytes_to_write);
        size_left -= bytes_to_write;
    }
}

inline char uchar_to_char(uint8_t c) {
    return *reinterpret_cast<char*>(&c);
}

inline uint8_t char_to_uchar(char c) {
    return *reinterpret_cast<uint8_t*>(&c);
}

template <typename pos_t>
inline char unsigned_to_char(pos_t c) {
    return *reinterpret_cast<char*>(&c);
}

template<typename T>
class no_init {
    static_assert(std::is_fundamental<T>::value);

    private:
    T v_;

    public:
    no_init() noexcept {}
    constexpr no_init(T value) noexcept: v_{value} {}
    constexpr operator T() const noexcept {return v_;}
};

template <typename T, typename Alloc = std::allocator<T>>
class default_init_allocator : public Alloc {
    using a_t = std::allocator_traits<Alloc>;
    using Alloc::Alloc;

    public:
    template<typename U> struct rebind {};
    template<typename U> void construct (U* ptr) noexcept(std::is_nothrow_default_constructible<U>::value) {::new(static_cast<void*>(ptr)) U;}
    template<typename U, typename... Args> void construct (U* ptr, Args&&... args) {}
};

void no_init_resize(std::string& str, size_t size) {
    (*reinterpret_cast<std::basic_string<char,std::char_traits<char>,default_init_allocator<char>>*>(&str)).resize(size);
}

template <typename T>
void no_init_resize(std::vector<T>& vec, size_t size) {
    (*reinterpret_cast<std::vector<no_init<T>>*>(&vec)).resize(size);
}

template <typename T1, typename T2>
void no_init_resize(std::vector<std::pair<T1,T2>>& vec, size_t size) {
    (*reinterpret_cast<std::vector<std::pair<no_init<T1>,no_init<T2>>>*>(&vec)).resize(size);
}

template <typename T>
void no_init_resize(std::vector<std::tuple<T,T,T>>& vec, size_t size) {
    (*reinterpret_cast<std::vector<std::tuple<no_init<T>,no_init<T>,no_init<T>>>*>(&vec)).resize(size);
}

void no_init_resize_with_exess(std::string& str, size_t size, size_t excess) {
    str.reserve(size + excess);
    no_init_resize(str, size);
    str.resize(size + excess);
    std::fill(str.end() - excess, str.end(), 0);
    str.resize(size);
}

template <int64_t start, int64_t end, int64_t inc, class T>
constexpr void for_constexpr(T&& f) {
    static_assert(inc != 0);
    if constexpr (inc > 0 ? start < end : start > end) {
        f(std::integral_constant<int64_t, start>());
        for_constexpr<start + inc, end, inc>(f);
    }
}

template <uint64_t start, uint64_t end, class T>
constexpr void for_constexpr_pow(T&& f) {
    static_assert(std::has_single_bit(start));
    static_assert(std::has_single_bit(end));
    static_assert(start <= end);
    f(std::integral_constant<uint64_t, start>());
    if constexpr (start < end) {
        for_constexpr_pow<2 * start, end>(f);
    }
}

template<class... Ts>
struct type_arr {
    template<std::uint8_t I>
    using get = std::tuple_element_t<I,std::tuple<Ts...>>;

    static constexpr std::uint8_t size = sizeof...(Ts);
};

template<uint8_t i, typename... Ts> using ith_type = typename std::tuple_element_t<i,std::tuple<Ts...>>;

template <uint8_t i, typename... Ts>
constexpr decltype(auto) ith_val(Ts&&... vals) noexcept {
    return std::get<i>(std::forward_as_tuple(vals...));
}

template <typename val_t, typename pos_t>
pos_t bin_search_min_geq(val_t value, pos_t left, pos_t right, std::function<val_t(pos_t)> value_at) {
    pos_t middle;

    while (left != right) {
        middle = left+(right-left)/2;

        if (value <= value_at(middle)) {
            right = middle;
        } else {
            left = middle+1;
        }
    }

    return left;
}

template <typename val_t, typename pos_t>
pos_t bin_search_max_lt(val_t value, pos_t left, pos_t right, std::function<val_t(pos_t)> value_at) {
    pos_t middle;

    while (left != right) {
        middle = left+(right-left)/2+1;

        if (value_at(middle) < value) {
            left = middle;
        } else {
            right = middle-1;
        }
    }

    return left;
}

template <typename val_t, typename pos_t>
pos_t bin_search_max_geq(val_t value, pos_t left, pos_t right, std::function<val_t(pos_t)> value_at) {
    pos_t middle;

    while (left != right) {
        middle = left+(right-left)/2+1;

        if (value_at(middle) >= value) {
            left = middle;
        } else {
            right = middle-1;
        }
    }

    return left;
}

template <typename val_t, typename pos_t>
pos_t bin_search_max_leq(val_t value, pos_t left, pos_t right, std::function<val_t(pos_t)> value_at) {
    pos_t middle;

    while (left != right) {
        middle = left+(right-left)/2+1;

        if (value_at(middle) <= value) {
            left = middle;
        } else {
            right = middle-1;
        }
    }

    return left;
}

template<typename int_t, typename... Ts>
constexpr int_t constexpr_sum(Ts&&... vals) {
    int_t sum = 0;

    for_constexpr<0, sizeof...(vals), 1>([&](auto i) {
        sum += ith_val<i>(vals...);
    });

    return sum;
}

std::string random_alphanumeric_string(uint64_t length) {
    static std::string possible_chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

    std::string str_rand;
    str_rand.reserve(length);
    
    for (uint64_t i=0; i<length; i++) {
        str_rand.push_back(possible_chars[std::rand()%possible_chars.size()]);
    }

    return str_rand;
}

enum direction {LEFT, RIGHT};

template <typename val_t, typename pos_t, direction search_dir>
pos_t exp_search_max_geq(val_t value, pos_t left, pos_t right, std::function<val_t(pos_t)> value_at) {
    if (right == left) {
        return left;
    }

    pos_t cur_step_size = 1;

    if constexpr (search_dir == LEFT) {
        right -= cur_step_size;

        while (value_at(right) < value) {
            cur_step_size *= 2;

            if (right < left + cur_step_size) {
                cur_step_size = right - left;
                right = left;
                break;
            }

            right -= cur_step_size;
        }

        return bin_search_max_geq<val_t, pos_t>(value, right, right + cur_step_size - 1, value_at);
    } else {
        left += cur_step_size;

        while (value_at(left) >= value) {
            cur_step_size *= 2;

            if (right < cur_step_size || right - cur_step_size < left) {
                cur_step_size = right - left + 1;
                left = right + 1;
                break;
            }

            left += cur_step_size;
        }

        return bin_search_max_geq<val_t, pos_t>(value, left - cur_step_size, left - 1, value_at);
    }
}

std::string format_query_throughput(uint64_t num_queries, uint64_t ns) {
    std::string str;
    double queries_per_ns = num_queries/(double)ns;

    if (queries_per_ns < 0.000001) {
        str = std::to_string(queries_per_ns*1000000000) + " queries/s";
    } else if (queries_per_ns < 0.001) {
        str = std::to_string(queries_per_ns*1000000) + " queries/ms";
    } else if (queries_per_ns < 1) {
        str = std::to_string(queries_per_ns*1000) + " queries/us";
    } else {
        str = std::to_string(queries_per_ns) + " queries/ns";
    }

    return str;
}

template<bool B, typename T>
struct constexpr_case {
    static constexpr bool value = B;
    using type = T;
};

template <bool B, typename TrueF, typename FalseF>
struct eval_if {
    using type = typename TrueF::type;
};

template <typename TrueF, typename FalseF>
struct eval_if<false, TrueF, FalseF> {
    using type = typename FalseF::type;
};    

template <bool B, typename T, typename F>
using eval_if_t = typename eval_if<B, T, F>::type;

template<typename Head, typename... Tail>
struct constexpr_switch {
    using type = eval_if_t<Head::value, Head, constexpr_switch<Tail...>>;
};

template <typename T>
struct constexpr_switch<T> {
    using type = T;
};

template <bool B, typename T>
struct constexpr_switch<constexpr_case<B, T>> {
    static_assert(B, "!");
    using type = T;
};

template<typename Head, typename... Tail>
using constexpr_switch_t = typename constexpr_switch<Head, Tail...>::type;

template<const auto &T, class = std::make_index_sequence<std::tuple_size<std::decay_t<decltype(T)>>::value>>
struct make_integer_sequence;

template<const auto &T, std::size_t ...Ts>
struct make_integer_sequence<T, std::index_sequence<Ts...>> {
    using type = std::integer_sequence<typename std::decay_t<decltype(T)>::value_type, T[Ts]...>;
};

template <typename int_t, int_t... seq>
constexpr int_t int_seq_at(std::integer_sequence<int_t, seq...>, std::size_t i) {
    return std::array<int_t, sizeof...(seq)>{seq...}[i];
}

template <typename uint_t>
inline static uint_t log2_clz(const uint_t val) {
    static_assert(std::is_same_v<uint_t, uint32_t> || std::is_same_v<uint_t, uint64_t>);
    
    if constexpr (std::is_same_v<uint_t, uint32_t>) {
        return 32 - __builtin_clz(val);
    } else {
        return 64 - __builtin_clzl(val);
    }
}

template <typename uint_t>
inline static uint_t log2_clz_up(const uint_t val) {
    static_assert(std::is_same_v<uint_t, uint32_t> || std::is_same_v<uint_t, uint64_t>);
    uint_t lz;

    if constexpr (std::is_same_v<uint_t, uint32_t>) {
        lz = __builtin_clz(val);
    } else {
        lz = __builtin_clzl(val);
    }

    return (sizeof(uint_t) * 8 - lz) + ((val << lz) != 0);
}

template <typename uint_t>
inline static uint_t div_ceil(const uint_t x, const uint_t y) {
    return x == 0 ? 0 : (1 + (x - 1) / y);
}

template <typename T>
T sum(std::vector<T> vec) {
    T sum = 0;
    for (T& v : vec) sum += v;
    return sum;
}

std::string random_repetitive_string(uint32_t min_size, uint32_t max_size) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
    std::uniform_int_distribution<char> char_distrib(
        std::numeric_limits<char>::min(),
        std::numeric_limits<char>::max());
    std::uniform_int_distribution<uint32_t> input_size_distrib(min_size, max_size);
    uint32_t target_input_size = input_size_distrib(mt);
    enum string_construction_operation {new_character = 0, repetition = 1, run = 2};
    double repetition_repetitiveness = prob_distrib(mt);
    double run_repetitiveness = prob_distrib(mt);
    std::uniform_int_distribution<uint32_t> repetition_length_distrib(
        1, (repetition_repetitiveness * target_input_size) / 100);
    std::uniform_int_distribution<uint32_t> run_length_distrib(
        1, (run_repetitiveness * target_input_size) / 200);
    std::discrete_distribution<uint8_t> next_operation_distrib({
        2 - (repetition_repetitiveness + run_repetitiveness),
        repetition_repetitiveness,
        run_repetitiveness
    });

    std::string input;
    no_init_resize_with_exess(input, target_input_size, 4 * 4096);
    input.clear();
    input.push_back(char_distrib(mt));

    while (input.size() < target_input_size) {
        switch (next_operation_distrib(mt)) {
            case new_character: input.push_back(char_distrib(mt)); break;
            case repetition: {
                uint32_t repetition_length = std::min<uint32_t>(
                    target_input_size - input.size(), repetition_length_distrib(mt));
                uint32_t repstition_source = std::rand()%input.size();
                for (uint32_t i = 0; i < repetition_length; i++)
                    input.push_back(input[repstition_source + i]);
                break;
            }
            case run: {
                uint32_t run_length = std::min<uint32_t>(
                    target_input_size - input.size(), run_length_distrib(mt));
                char run_char = char_distrib(mt);
                for (uint32_t i = 0; i < run_length; i++)
                    input.push_back(run_char);
                break;
            }
        }
    }

    return input;
}