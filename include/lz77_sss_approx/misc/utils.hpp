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

#include <malloc_count/malloc_count.h>

using uint128_t = alx::rolling_hash::uint128_t;

std::string random_repetitive_string(uint32_t max_size) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
    std::uniform_int_distribution<char> char_distrib(-128, 127);
    uint32_t target_input_size = std::max<uint32_t>(max_size, 10000);
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
    input.reserve(target_input_size);
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

std::string format_threads(uint16_t p) {
    if (p == 1) {
        return "1 thread";
    } else {
        return std::to_string(p) + " threads";
    }
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

template <int64_t start, int64_t end, int64_t inc, class T>
constexpr void for_constexpr(T&& f) {
    static_assert(inc != 0);

    if constexpr (inc > 0) {
        if constexpr (start < end) {
            f(std::integral_constant<int64_t, start>());
            for_constexpr<start + inc, end, inc>(f);
        }
    } else {
        if constexpr (start > end) {
            f(std::integral_constant<int64_t, start>());
            for_constexpr<start + inc, end, inc>(f);
        }
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

template <typename pos_t, typename val_t>
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