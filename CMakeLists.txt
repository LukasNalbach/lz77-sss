cmake_minimum_required(VERSION "3.10.2")
set (CMAKE_CXX_STANDARD 20)
project (lz77_sss LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(lz77_sss_debug_flags "-Wall -Wextra -Wpedantic -ggdb")
set(lz77_sss_release_flags "-w -march=native -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG           "${lz77_sss_debug_flags} -O0")
set(CMAKE_CXX_FLAGS_RELEASE         "${lz77_sss_release_flags} -Ofast")
set(CMAKE_CXX_FLAGS_MINSIZEREL      "${lz77_sss_release_flags} -Os")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${lz77_sss_debug_flags} -Ofast")

option(LZ77_SSS_BUILD_CLI "Build the cli programs" ON)
option(LZ77_SSS_BUILD_BENCH_CLI "Build the benchmark cli program" ON)
option(LZ77_SSS_BUILD_EXAMPLES "Build the example programs" ON)
option(LZ77_SSS_BUILD_TESTS "Build tests" ON)

add_subdirectory("${CMAKE_SOURCE_DIR}/external/lz77/")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/alx/")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/fp/")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/sparse-map/")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/gtl/")

set_target_properties(
  ips4o_example benchmark_lce benchmark_pred benchmark_ss
  gen_sss genqueries test_lce test_pred test_rmq
  test_rolling_hash test_string_synchronizing_set
  PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1
)

# lz77_sss
add_library(lz77_sss INTERFACE)
target_include_directories(lz77_sss INTERFACE "${CMAKE_SOURCE_DIR}/include/")
target_link_libraries(lz77_sss INTERFACE malloc_count alx_rolling_hash alx_lce lz77 fp sparse_map gtl)

if(LZ77_SSS_BUILD_CLI)
  add_executable(lz77_sss_3-aprx "${CMAKE_SOURCE_DIR}/cli/lz77_sss_3-aprx.cpp")
  target_link_libraries(lz77_sss_3-aprx PRIVATE lz77_sss)

  add_executable(lz77_sss_1_5-aprx "${CMAKE_SOURCE_DIR}/cli/lz77_sss_1_5-aprx.cpp")
  target_link_libraries(lz77_sss_1_5-aprx PRIVATE lz77_sss)

  add_executable(lz77_sss_exact "${CMAKE_SOURCE_DIR}/cli/lz77_sss_exact.cpp")
  target_link_libraries(lz77_sss_exact PRIVATE lz77_sss)

  add_executable(lz77_sss_exact-smpl "${CMAKE_SOURCE_DIR}/cli/lz77_sss_exact-smpl.cpp")
  target_link_libraries(lz77_sss_exact-smpl PRIVATE lz77_sss)

  add_executable(lz77_sss_decode "${CMAKE_SOURCE_DIR}/cli/lz77_sss_decode.cpp")
  target_link_libraries(lz77_sss_decode PRIVATE lz77_sss)
endif()

if(LZ77_SSS_BUILD_BENCH_CLI)
  add_executable(lz77_sss_bench "${CMAKE_SOURCE_DIR}/cli/lz77_sss_bench.cpp")
  target_link_libraries(lz77_sss_bench PRIVATE lz77_sss)

  add_executable(lz77_sss_bench-tau "${CMAKE_SOURCE_DIR}/cli/lz77_sss_bench-tau.cpp")
  target_link_libraries(lz77_sss_bench-tau PRIVATE lz77_sss)
endif()

if(LZ77_SSS_BUILD_EXAMPLES)
  add_executable(lz77_sss_example "${CMAKE_SOURCE_DIR}/examples/example.cpp")
  target_link_libraries(lz77_sss_example PRIVATE lz77_sss)
endif()

if(LZ77_SSS_BUILD_TESTS)
  enable_testing()
  include(GoogleTest)
  include(Dart)

  add_executable(test-lz77_sss "${CMAKE_SOURCE_DIR}/tests/test_lz77_sss.cpp")
  target_link_libraries(test-lz77_sss PRIVATE gtest_main lz77_sss)
  gtest_discover_tests(test-lz77_sss)

  add_executable(test-sample_index "${CMAKE_SOURCE_DIR}/tests/test_sample_index.cpp")
  target_link_libraries(test-sample_index PRIVATE gtest_main lz77_sss)
  gtest_discover_tests(test-sample_index)

  add_executable(test-dynamic_range "${CMAKE_SOURCE_DIR}/tests/test_dynamic_range.cpp")
  target_link_libraries(test-dynamic_range PRIVATE gtest_main lz77_sss)
  gtest_discover_tests(test-dynamic_range)

  add_executable(test-static_weighted_range "${CMAKE_SOURCE_DIR}/tests/test_static_weighted_range.cpp")
  target_link_libraries(test-static_weighted_range PRIVATE gtest_main lz77_sss)
  gtest_discover_tests(test-static_weighted_range)

  add_executable(test-decomposed_range "${CMAKE_SOURCE_DIR}/tests/test_decomposed_range.cpp")
  target_link_libraries(test-decomposed_range PRIVATE gtest_main lz77_sss)
  gtest_discover_tests(test-decomposed_range)

  add_executable(test-rk61_substring "${CMAKE_SOURCE_DIR}/tests/test_rk61_substring.cpp")
  target_link_libraries(test-rk61_substring PRIVATE gtest_main lz77_sss)
  gtest_discover_tests(test-rk61_substring)
endif()