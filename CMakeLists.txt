cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 23)
project (lz77_sss LANGUAGES C CXX)

option(LZ77_SSS_BUILD_CLI           "Build cli programs"            ON)
option(LZ77_SSS_BUILD_BENCH         "Build benchmark cli programs"  OFF)
option(LZ77_SSS_BUILD_TESTS         "Build tests"                   ON)
option(LZ77_SSS_EXTERN_TEMPLATES    "Compile the factorizer once instead of in every translation unit" ON)
set(LZ77_SSS_LCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/lce" CACHE PATH "lce source directory")

if(WIN32)
  set(LZ77_SSS_BUILD_BENCH OFF CACHE BOOL "" FORCE)
endif()

# add external/string-utils, if the parent project has not added it
if(NOT TARGET string_utils)
  file(GLOB lz77_sss_string_utils_contents "${CMAKE_CURRENT_SOURCE_DIR}/external/string-utils/*")

  if(lz77_sss_string_utils_contents STREQUAL "")
    find_package(Git REQUIRED)
    message(STATUS "Downloading the git submodule external/string-utils")
    execute_process(COMMAND "${GIT_EXECUTABLE}" submodule update --init -- external/string-utils
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()

  # the prefix of the options and environment variables (LZ77_SSS_MARCH_NATIVE, LZ77_SSS_TEST_ITERATIONS, ...)
  if(PROJECT_IS_TOP_LEVEL)
    set(STRING_UTILS_PREFIX LZ77_SSS)
    option(LZ77_SSS_ENABLE_LTO "Enable link-time optimization (LTO/IPO) for lz77-sss targets" ON)
    include("${CMAKE_CURRENT_SOURCE_DIR}/external/string-utils/cmake/settings.cmake")
    set(CMAKE_CXX_STANDARD 23)
  endif()

  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/string-utils" EXCLUDE_FROM_ALL)
endif()

# the libraries of string-utils beyond its core that are used here
string_utils_add_components(libsais)

# the text classes
if(NOT TARGET compact_texts)
  string_utils_fetch_submodules("${CMAKE_CURRENT_SOURCE_DIR}" PATHS external/compact-texts)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/compact-texts" EXCLUDE_FROM_ALL)
endif()

if(NOT EXISTS "${LZ77_SSS_LCE_DIR}/CMakeLists.txt")
  string_utils_fetch_submodules("${CMAKE_CURRENT_SOURCE_DIR}" PATHS external/lce)
endif()

set(CMAKE_POLICY_VERSION_MINIMUM 3.10)
add_subdirectory("${LZ77_SSS_LCE_DIR}" "${CMAKE_BINARY_DIR}/external/lce" EXCLUDE_FROM_ALL)
unset(CMAKE_POLICY_VERSION_MINIMUM)

function(lz77_sss_apply_patch patch dir)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${patch}")
  get_filename_component(_parent "${dir}" DIRECTORY)
  set(_git_apply "${CMAKE_COMMAND}" -E env "GIT_CEILING_DIRECTORIES=${_parent}" "${GIT_EXECUTABLE}" apply)
  execute_process(COMMAND ${_git_apply} --reverse --check "${patch}"
    WORKING_DIRECTORY "${dir}" RESULT_VARIABLE _applied OUTPUT_QUIET ERROR_QUIET)
  if(NOT _applied EQUAL 0)
    execute_process(COMMAND ${_git_apply} "${patch}"
      WORKING_DIRECTORY "${dir}" RESULT_VARIABLE _result ERROR_VARIABLE _error)
    if(NOT _result EQUAL 0)
      message(FATAL_ERROR "${patch} does not apply to ${dir}:\n${_error}")
    endif()
    message(STATUS "Applied ${patch}")
  endif()
endfunction()

lz77_sss_apply_patch("${CMAKE_CURRENT_SOURCE_DIR}/patches/lce.patch" "${LZ77_SSS_LCE_DIR}")

add_library(lz77_sss_core INTERFACE)
target_include_directories(lz77_sss_core INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_link_libraries(lz77_sss_core INTERFACE string_utils libsais libsais40 libsais48 compact_texts lce)

add_library(lz77_sss INTERFACE)
target_link_libraries(lz77_sss INTERFACE lz77_sss_core)

if(LZ77_SSS_EXTERN_TEMPLATES)
  add_library(lz77_sss_inst STATIC)
  target_link_libraries(lz77_sss_inst PUBLIC lz77_sss_core)

  foreach(text direct_text packed_text split_text int_direct_text int64_direct_text int_packed_text int_split_text)
    foreach(part factorizer lce lpf fact sa aprx_to_exact)
      set(obj lz77_sss_${text}_${part})
      add_library(${obj} OBJECT "${CMAKE_CURRENT_SOURCE_DIR}/src/inst/${part}.cpp")
      target_link_libraries(${obj} PRIVATE lz77_sss_core)
      target_compile_definitions(${obj} PRIVATE LZ77_SSS_INST_TEXT=${text})
      target_compile_options(${obj} PRIVATE $<$<CXX_COMPILER_ID:GNU>:--param=inline-unit-growth=200>)
      target_sources(lz77_sss_inst PRIVATE $<TARGET_OBJECTS:${obj}>)
    endforeach()
  endforeach()

  target_link_libraries(lz77_sss INTERFACE lz77_sss_inst)
else()
  target_compile_definitions(lz77_sss INTERFACE LZ77_SSS_NO_EXTERN_TEMPLATES)
endif()

function(lz77_sss_add_executable name src outdir)
  string_utils_add_executable(${name} "${src}" "${outdir}" lz77_sss ${ARGN})
endfunction()

if(LZ77_SSS_BUILD_CLI)
  foreach(tool ssszip lz77-sss-3-aprx lz77-sss-exact
               lz77-sss-exact-smpl lz77-sss-decode)
    string(REPLACE "-" "_" src "${tool}")
    lz77_sss_add_executable(${tool} "${CMAKE_CURRENT_SOURCE_DIR}/cli/${src}.cpp" "${CMAKE_BINARY_DIR}/cli/")
  endforeach()
endif()

if(LZ77_SSS_BUILD_BENCH)
  foreach(tool lz77-sss-bench lz77-sss-bench-tau gen-range-queries bench-range-queries zip-bench)
    string(REPLACE "-" "_" src "${tool}")
    lz77_sss_add_executable(${tool} "${CMAKE_CURRENT_SOURCE_DIR}/bench/${src}.cpp" "${CMAKE_BINARY_DIR}/bench/")
  endforeach()

  string_utils_fetch_submodules("${CMAKE_CURRENT_SOURCE_DIR}" PATHS external/alz)
  include(ExternalProject)
  set(lz77_sss_alz_binary "${CMAKE_BINARY_DIR}/bench/alz")

  ExternalProject_Add(alz_cli
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/alz"
    PREFIX "${CMAKE_BINARY_DIR}/alz"
    CMAKE_ARGS
      "-DCMAKE_BUILD_TYPE=Release"
      "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
      "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target alz-cli
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "<BINARY_DIR>/src/alz" "${lz77_sss_alz_binary}"
    BUILD_BYPRODUCTS "${lz77_sss_alz_binary}"
  )

  add_dependencies(zip-bench alz_cli)
endif()

if(LZ77_SSS_BUILD_TESTS)
  include(GoogleTest)
  enable_testing()
  string_utils_add_components(tests)

  foreach(s lz77-sss sample-index dynamic-range static-weighted-range
            decomposed-range rabin-karp-substring)
    string(REPLACE "-" "_" src "test-${s}")
    lz77_sss_add_executable("test-${s}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/${src}.cpp"
      "${CMAKE_BINARY_DIR}/tests/" string_utils_test)
    gtest_discover_tests("test-${s}")
  endforeach()
endif()
