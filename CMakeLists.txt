cmake_minimum_required(VERSION "3.10.2")
set (CMAKE_CXX_STANDARD 23)
project (lz77_sss LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message("Building in ${CMAKE_BUILD_TYPE} mode")

option(LZ77_SSS_BUILD_CLI           "Build cli programs"            ON)
option(LZ77_SSS_BUILD_BENCH         "Build benchmark cli programs"  ON)
option(LZ77_SSS_BUILD_EXAMPLES      "Build example programs"        OFF)
option(LZ77_SSS_BUILD_TESTS         "Build tests"                   ON)
option(LZ77_SSS_USE_MALLOC_COUNT    "Track memory usage via malloc_count (overrides malloc/free)" ON)
option(LZ77_SSS_MARCH_NATIVE        "Tune the build for the host CPU (-march=native); turn OFF for portable binaries" ON)
option(LZ77_SSS_ENABLE_LTO          "Enable link-time optimization (LTO/IPO) for lz77-sss targets" ON)

if(WIN32)
  set(LZ77_SSS_BUILD_BENCH OFF CACHE BOOL "" FORCE)
endif()

if(LZ77_SSS_MARCH_NATIVE)
  set(lz77_sss_arch_flags "-march=native")
else()
  set(lz77_sss_arch_flags "")
endif()

add_compile_options(-pthread)

set(lz77_sss_debug_flags            "-Wall -Wextra -Wpedantic -ggdb")
set(lz77_sss_release_flags          "-w ${lz77_sss_arch_flags} -DNDEBUG")

set(CMAKE_CXX_FLAGS_DEBUG           "${lz77_sss_debug_flags}        -O0")
set(CMAKE_CXX_FLAGS_RELEASE         "${lz77_sss_release_flags}      -O3")
set(CMAKE_CXX_FLAGS_MINSIZEREL      "${lz77_sss_release_flags}      -Os")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${lz77_sss_debug_flags}        -O3")

if(WIN32)
  include_directories(BEFORE "${CMAKE_SOURCE_DIR}/cmake/winshim")
  add_compile_definitions(_REENTRANT)
  set(IPS4O_OPTIMIZE_FOR_NATIVE OFF CACHE BOOL "" FORCE)
  set(IPS4O_DISABLE_PARALLEL ON CACHE BOOL "" FORCE)

  if(NOT MINGW)
    add_compile_definitions(MSVC_COMPILER _HAS_DEPRECATED_RESULT_OF=1)
  else()
    get_filename_component(_mingw_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
    get_filename_component(_mingw_root "${_mingw_bin}" DIRECTORY)
    if(EXISTS "${_mingw_root}/lib/libgomp.a")
      set(OpenMP_gomp_LIBRARY "${_mingw_root}/lib/libgomp.a" CACHE FILEPATH "" FORCE)
    endif()
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-mcx16)
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_definitions(HURCHALLA_DISALLOW_INLINE_ASM_MODMUL)
    add_link_options("LINKER:/defaultlib:libcmt,/defaultlib:libcpmt,/defaultlib:libvcruntime,/defaultlib:libucrt,/defaultlib:oldnames")
    add_link_options("LINKER:/STACK:0x10000000")
    get_filename_component(_llvm_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
    get_filename_component(_llvm_root "${_llvm_bin}" DIRECTORY)
    file(GLOB _lz77_clang_rt "${_llvm_root}/lib/clang/*/lib/windows/clang_rt.builtins-x86_64.lib")
    if(_lz77_clang_rt)
      list(GET _lz77_clang_rt 0 LZ77_SSS_WIN_CLANG_RT)
    endif()
    if(EXISTS "${_llvm_root}/lib/libomp.lib")
      foreach(_lang C CXX)
        set(OpenMP_${_lang}_FLAGS "-fopenmp=libomp" CACHE STRING "" FORCE)
        set(OpenMP_${_lang}_LIB_NAMES "libomp" CACHE STRING "" FORCE)
      endforeach()
      set(OpenMP_libomp_LIBRARY "${_llvm_root}/lib/libomp.lib" CACHE FILEPATH "" FORCE)
      if(EXISTS "${_llvm_root}/bin/libomp.dll")
        list(APPEND LZ77_SSS_WIN_RUNTIME_DLLS "${_llvm_root}/bin/libomp.dll")
      endif()
    endif()
  endif()

  set(LZ77_SSS_WINDOWS_TBB_ROOT "" CACHE PATH "Root of an extracted oneAPI TBB (Windows) release")
  if(LZ77_SSS_WINDOWS_TBB_ROOT)
    include_directories(BEFORE "${LZ77_SSS_WINDOWS_TBB_ROOT}/include")
    find_library(LZ77_SSS_WIN_TBB_LIB NAMES tbb12 tbb
      PATHS "${LZ77_SSS_WINDOWS_TBB_ROOT}/lib/intel64/vc14" "${LZ77_SSS_WINDOWS_TBB_ROOT}/lib" NO_DEFAULT_PATH)
    file(GLOB _tbb_dll "${LZ77_SSS_WINDOWS_TBB_ROOT}/redist/intel64/vc14/tbb12.dll"
                       "${LZ77_SSS_WINDOWS_TBB_ROOT}/redist/intel64/vc14/tbb.dll"
                       "${LZ77_SSS_WINDOWS_TBB_ROOT}/bin/*tbb12*.dll"
                       "${LZ77_SSS_WINDOWS_TBB_ROOT}/bin/*tbb.dll")
    if(_tbb_dll)
      list(APPEND LZ77_SSS_WIN_RUNTIME_DLLS ${_tbb_dll})
    endif()
  else()
    message(WARNING "LZ77_SSS_WINDOWS_TBB_ROOT is not set: ips4o's parallel sort needs Intel "
      "oneTBB on Windows. Download an 'oneapi-tbb-*-win' release and set LZ77_SSS_WINDOWS_TBB_ROOT.")
  endif()
endif()

find_package(OpenMP REQUIRED)

set(LCE_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(LCE_USE_MALLOC_COUNT ${LZ77_SSS_USE_MALLOC_COUNT} CACHE BOOL "" FORCE)
set(LCE_MARCH_NATIVE ${LZ77_SSS_MARCH_NATIVE} CACHE BOOL "" FORCE)
if(WIN32 AND LZ77_SSS_WINDOWS_TBB_ROOT)
  set(LCE_WINDOWS_TBB_ROOT "${LZ77_SSS_WINDOWS_TBB_ROOT}" CACHE PATH "" FORCE)
endif()

set(CMAKE_POLICY_VERSION_MINIMUM 3.10)
add_subdirectory("${CMAKE_SOURCE_DIR}/external/lce/")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/fp/")
add_subdirectory("${CMAKE_SOURCE_DIR}/external/sparse-map/")
unset(CMAKE_POLICY_VERSION_MINIMUM)

set_target_properties(ips4o_example PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

foreach(_la_target example tests)
  if(TARGET ${_la_target})
    set_target_properties(${_la_target} PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
  endif()
endforeach()

add_library(lz77_sss INTERFACE)
target_include_directories(lz77_sss INTERFACE "${CMAKE_SOURCE_DIR}/include/")
target_link_libraries(lz77_sss INTERFACE malloc_count lce libsais fp sparse_map OpenMP::OpenMP_CXX)

if(WIN32 AND LZ77_SSS_WIN_TBB_LIB)
  target_link_libraries(lz77_sss INTERFACE "${LZ77_SSS_WIN_TBB_LIB}")
endif()

if(WIN32 AND LZ77_SSS_WIN_CLANG_RT)
  target_link_libraries(lz77_sss INTERFACE "${LZ77_SSS_WIN_CLANG_RT}")
endif()

if(MINGW)
  target_link_libraries(lz77_sss INTERFACE atomic dl)
  target_link_options(lz77_sss INTERFACE -static)
endif()

if(LZ77_SSS_ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT LZ77_SSS_IPO_SUPPORTED OUTPUT LZ77_SSS_IPO_ERROR LANGUAGES CXX)
  if(NOT LZ77_SSS_IPO_SUPPORTED)
    message(WARNING "LZ77_SSS_ENABLE_LTO is ON but LTO/IPO is not supported by the toolchain: ${LZ77_SSS_IPO_ERROR}")
  endif()
endif()

function(lz77_sss_add_tool name src outdir)
  add_executable(${name} "${src}")
  target_link_libraries(${name} PRIVATE lz77_sss ${ARGN})
  set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${outdir}")
  if(LZ77_SSS_ENABLE_LTO AND LZ77_SSS_IPO_SUPPORTED)
    set_target_properties(${name} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
  endif()
  foreach(_dll ${LZ77_SSS_WIN_RUNTIME_DLLS})
    add_custom_command(TARGET ${name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "${_dll}" "${outdir}")
  endforeach()
endfunction()

if(LZ77_SSS_BUILD_CLI)
  foreach(tool ssszip lz77-sss-3-aprx lz77-sss-lpf-lnf-aprx lz77-sss-exact
               lz77-sss-exact-smpl lz77-sss-decode)
    string(REPLACE "-" "_" src "${tool}")
    lz77_sss_add_tool(${tool} "${CMAKE_SOURCE_DIR}/cli/${src}.cpp" "${CMAKE_BINARY_DIR}/cli/")
  endforeach()

  if(EXISTS "${CMAKE_SOURCE_DIR}/cli/lz77_sss_ref_analysis.cpp")
    lz77_sss_add_tool(lz77-sss-ref-analysis "${CMAKE_SOURCE_DIR}/cli/lz77_sss_ref_analysis.cpp"
      "${CMAKE_BINARY_DIR}/cli/")
  endif()
endif()

if(LZ77_SSS_BUILD_BENCH)
  foreach(tool lz77-sss-bench lz77-sss-bench-tau gen-range-queries bench-range-queries zip-bench)
    string(REPLACE "-" "_" src "${tool}")
    lz77_sss_add_tool(${tool} "${CMAKE_SOURCE_DIR}/bench/${src}.cpp" "${CMAKE_BINARY_DIR}/bench/")
  endforeach()
endif()

if(LZ77_SSS_BUILD_EXAMPLES)
  lz77_sss_add_tool(lz77-sss-example "${CMAKE_SOURCE_DIR}/examples/lz77_sss_example.cpp"
    "${CMAKE_BINARY_DIR}/examples/")
endif()

if(LZ77_SSS_BUILD_TESTS)
  if(NOT TARGET gtest)
    add_subdirectory("${CMAKE_SOURCE_DIR}/external/googletest/")
  endif()
  include(GoogleTest)
  enable_testing()

  foreach(_gtest_target gtest gtest_main gmock gmock_main)
    if(TARGET ${_gtest_target})
      set_target_properties(${_gtest_target}
        PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "")
    endif()
  endforeach()

  foreach(s lz77-sss sample-index dynamic-range static-weighted-range
            decomposed-range rabin-karp-substring)
    string(REPLACE "-" "_" src "test-${s}")
    lz77_sss_add_tool("test-${s}" "${CMAKE_SOURCE_DIR}/tests/${src}.cpp"
      "${CMAKE_BINARY_DIR}/tests/" gtest_main)
    gtest_discover_tests("test-${s}")
  endforeach()
endif()
