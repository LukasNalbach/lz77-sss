cmake_minimum_required(VERSION "3.10.2")
set (CMAKE_CXX_STANDARD 20)
project (lz77_sss_approx LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(lz77_sss_approx_debug_flags "-Wall -Wextra -Wpedantic -ggdb")
set(lz77_sss_approx_release_flags "-w -march=native -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG           "${lz77_sss_approx_debug_flags} -O0")
set(CMAKE_CXX_FLAGS_RELEASE         "${lz77_sss_approx_release_flags} -Ofast")
set(CMAKE_CXX_FLAGS_MINSIZEREL      "${lz77_sss_approx_release_flags} -Os")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${lz77_sss_approx_debug_flags} -Ofast")

option(LZ77_SSS_APRX_BUILD_BENCH_CLI "Build the benchmark cli program" ON)
option(LZ77_SSS_APRX_BUILD_EXAMPLES "Build the example programs" ON)
option(LZ77_SSS_APRX_BUILD_TESTS "Build tests" ON)

# abseil
set(ABSL_PROPAGATE_CXX_STD ON CACHE BOOL "" FORCE)
add_subdirectory("${CMAKE_SOURCE_DIR}/external/abseil-cpp/")

# lz77
add_subdirectory("${CMAKE_SOURCE_DIR}/external/lz77/")

# alx
add_subdirectory("${CMAKE_SOURCE_DIR}/external/alx/")

# fp
add_subdirectory("${CMAKE_SOURCE_DIR}/external/fp/")

set_target_properties(
  ips4o_example benchmark_lce benchmark_pred benchmark_ss
  gen_sss genqueries test_lce test_pred test_rmq
  test_rolling_hash test_string_synchronizing_set
  PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1
)

# lz77_sss_approx
add_library(lz77_sss_approx INTERFACE)
target_include_directories(lz77_sss_approx INTERFACE "${CMAKE_SOURCE_DIR}/include/")
target_link_libraries(lz77_sss_approx INTERFACE ips4o malloc_count absl::btree alx_rolling_hash alx_lce lz77 fp)
# semidynamic-compact-index OnlineRlbwt gtl

if(LZ77_SSS_APRX_BUILD_BENCH_CLI)
  add_executable(lz77_sss_approx_bench "${CMAKE_SOURCE_DIR}/cli/bench.cpp")
  target_link_libraries(lz77_sss_approx_bench PRIVATE lz77_sss_approx)
endif()

if(LZ77_SSS_APRX_BUILD_EXAMPLES)
  add_executable(lz77_sss_approx_example "${CMAKE_SOURCE_DIR}/examples/example.cpp")
  target_link_libraries(lz77_sss_approx_example PRIVATE lz77_sss_approx)
endif()

if(LZ77_SSS_APRX_BUILD_TESTS)
  enable_testing()
  include(GoogleTest)
  include(Dart)

  add_executable(test-lz77_sss_aprx "${CMAKE_SOURCE_DIR}/tests/test.cpp")
  target_link_libraries(test-lz77_sss_aprx PRIVATE gtest_main lz77_sss_approx)
  gtest_discover_tests(test-lz77_sss_aprx)
endif()