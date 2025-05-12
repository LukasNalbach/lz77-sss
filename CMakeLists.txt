cmake_minimum_required(VERSION "3.10.2")
set (CMAKE_CXX_STANDARD 20)
project (lz77_sss LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(lz77_sss_debug_flags            "-Wall -Wextra -Wpedantic -ggdb")
set(lz77_sss_release_flags          "-w -march=native -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG           "${lz77_sss_debug_flags}        -O0")
set(CMAKE_CXX_FLAGS_RELEASE         "${lz77_sss_release_flags}      -Ofast")
set(CMAKE_CXX_FLAGS_MINSIZEREL      "${lz77_sss_release_flags}      -Os")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${lz77_sss_debug_flags}        -Ofast")

option(LZ77_SSS_BUILD_CLI           "Build cli programs"            ON)
option(LZ77_SSS_BUILD_BENCH         "Build benchmark cli programs"  OFF)
option(LZ77_SSS_BUILD_EXAMPLES      "Build example programs"        OFF)
option(LZ77_SSS_BUILD_TESTS         "Build tests"                   OFF)

find_package(OpenMP REQUIRED)

add_subdirectory("${CMAKE_SOURCE_DIR}/external/lz77/")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/lce/")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/fp/")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/sparse-map/")

set_target_properties(
  ips4o_example benchmark_lce benchmark_pred benchmark_lce
  gen_sss gen_queries gen_sa_lcp test_lce test_pred test_rmq
  test_rolling_hash test_string_synchronizing_set
  PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

add_library(lz77_sss INTERFACE)
target_include_directories(lz77_sss INTERFACE "${CMAKE_SOURCE_DIR}/include/")
target_link_libraries(lz77_sss INTERFACE malloc_count lce lz77 fp sparse_map OpenMP::OpenMP_CXX)

if(LZ77_SSS_BUILD_CLI)
  add_executable(ssszip "${CMAKE_SOURCE_DIR}/cli/ssszip.cpp")
  target_link_libraries(ssszip PRIVATE lz77_sss)

  add_executable(lz77-sss-3-aprx "${CMAKE_SOURCE_DIR}/cli/lz77_sss_3_aprx.cpp")
  target_link_libraries(lz77-sss-3-aprx PRIVATE lz77_sss)

  add_executable(lz77-sss-lpf-lnf-aprx "${CMAKE_SOURCE_DIR}/cli/lz77_sss_lpf_lnf_aprx.cpp")
  target_link_libraries(lz77-sss-lpf-lnf-aprx PRIVATE lz77_sss)

  add_executable(lz77-sss-exact "${CMAKE_SOURCE_DIR}/cli/lz77_sss_exact.cpp")
  target_link_libraries(lz77-sss-exact PRIVATE lz77_sss)

  add_executable(lz77-sss-exact-smpl "${CMAKE_SOURCE_DIR}/cli/lz77_sss_exact_smpl.cpp")
  target_link_libraries(lz77-sss-exact-smpl PRIVATE lz77_sss)

  add_executable(lz77-sss-decode "${CMAKE_SOURCE_DIR}/cli/lz77_sss_decode.cpp")
  target_link_libraries(lz77-sss-decode PRIVATE lz77_sss)

  set_target_properties(ssszip lz77-sss-3-aprx lz77-sss-lpf-lnf-aprx
    lz77-sss-exact lz77-sss-exact-smpl lz77-sss-decode
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/cli/")
endif()

if(LZ77_SSS_BUILD_BENCH)
  add_executable(lz77-sss-bench "${CMAKE_SOURCE_DIR}/bench/lz77_sss_bench.cpp")
  target_link_libraries(lz77-sss-bench PRIVATE lz77_sss)

  add_executable(lz77-sss-bench-tau "${CMAKE_SOURCE_DIR}/bench/lz77_sss_bench_tau.cpp")
  target_link_libraries(lz77-sss-bench-tau PRIVATE lz77_sss)

  add_executable(gen-range-queries "${CMAKE_SOURCE_DIR}/bench/gen_range_queries.cpp")
  target_link_libraries(gen-range-queries PRIVATE lz77_sss)

  add_executable(bench-range-queries "${CMAKE_SOURCE_DIR}/bench/bench_range_queries.cpp")
  target_link_libraries(bench-range-queries PRIVATE lz77_sss)

  add_executable(zip-bench "${CMAKE_SOURCE_DIR}/bench/zip_bench.cpp")
  target_link_libraries(zip-bench PRIVATE lz77_sss)

  set_target_properties(lz77-sss-bench lz77-sss-bench-tau
    gen-range-queries bench-range-queries zip-bench
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bench/")
endif()

if(LZ77_SSS_BUILD_EXAMPLES)
  add_executable(lz77-sss-example "${CMAKE_SOURCE_DIR}/examples/lz77_sss_example.cpp")
  target_link_libraries(lz77-sss-example PRIVATE lz77_sss)
  set_target_properties(lz77-sss-example PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/examples/")
endif()

if(LZ77_SSS_BUILD_TESTS)
  enable_testing()
  find_package(GTest REQUIRED)
  include_directories(${GTEST_INCLUDE_DIR})

  add_executable(test-lz77-sss "${CMAKE_SOURCE_DIR}/tests/test_lz77_sss.cpp")
  target_link_libraries(test-lz77-sss PRIVATE gtest_main lz77_sss)

  add_executable(test-sample-index "${CMAKE_SOURCE_DIR}/tests/test_sample_index.cpp")
  target_link_libraries(test-sample-index PRIVATE gtest_main lz77_sss)

  add_executable(test-dynamic-range "${CMAKE_SOURCE_DIR}/tests/test_dynamic_range.cpp")
  target_link_libraries(test-dynamic-range PRIVATE gtest_main lz77_sss)

  add_executable(test-static-weighted-range "${CMAKE_SOURCE_DIR}/tests/test_static_weighted_range.cpp")
  target_link_libraries(test-static-weighted-range PRIVATE gtest_main lz77_sss)

  add_executable(test-decomposed-range "${CMAKE_SOURCE_DIR}/tests/test_decomposed_range.cpp")
  target_link_libraries(test-decomposed-range PRIVATE gtest_main lz77_sss)

  add_executable(test-rabin-karp-substring "${CMAKE_SOURCE_DIR}/tests/test_rabin_karp_substring.cpp")
  target_link_libraries(test-rabin-karp-substring PRIVATE gtest_main lz77_sss)

  gtest_discover_tests(test-lz77-sss test-sample-index test-dynamic-range
    test-static-weighted-range test-decomposed-range test-rabin-karp-substring)

  set_target_properties(test-lz77-sss test-sample-index test-dynamic-range
    test-static-weighted-range test-decomposed-range test-rabin-karp-substring
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
endif()