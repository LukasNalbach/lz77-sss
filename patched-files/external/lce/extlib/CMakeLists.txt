set(ALX_EXTLIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS "alx extlib path: ${ALX_EXTLIB_DIR}")

add_library(rk-lce INTERFACE)
target_include_directories(rk-lce INTERFACE rk-lce/internal) 

add_library(tlx_clp tlx/tlx/cmdline_parser.cpp tlx/tlx/string/parse_si_iec_units.cpp)
target_include_directories(tlx_clp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tlx)

add_subdirectory(fmt)

if(LZ77_SSS_MALLOC_COUNT_STUB AND (NOT LZ77_SSS_USE_MALLOC_COUNT OR WIN32))
  add_library(malloc_count "${LZ77_SSS_MALLOC_COUNT_STUB}")
  target_compile_definitions(malloc_count PUBLIC LZ77_SSS_USE_MALLOC_COUNT=0)
else()
  add_library(malloc_count malloc_count/malloc_count.c)
  target_link_libraries(malloc_count ${CMAKE_DL_LIBS})
  target_compile_definitions(malloc_count PUBLIC LZ77_SSS_USE_MALLOC_COUNT=1)
endif()
target_include_directories(malloc_count PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/modular_arithmetic)

add_library(parallel-hashmap INTERFACE)
target_include_directories(parallel-hashmap INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/parallel-hashmap/)

add_library(pgm INTERFACE)
target_link_libraries(pgm INTERFACE  OpenMP::OpenMP_CXX)
target_include_directories(pgm INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/PGM-index/include/pgm)

# ips4o/extern/cmake-modules messes up loading TBB
option(KEEP_IPS4O_EXTERN_FOLDER OFF)
if(NOT ${KEEP_IPS4O_EXTERN_FOLDER})
file(REMOVE_RECURSE ${CMAKE_CURRENT_SOURCE_DIR}/ips4o/extern/cmake-modules)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ips4o/extern/cmake-modules)
endif()
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ips4o)

add_library(gsaca_ds INTERFACE)
target_link_libraries(gsaca_ds INTERFACE OpenMP::OpenMP_CXX ips4o -latomic)
target_include_directories(gsaca_ds INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/gsaca-double-sort)

add_subdirectory(la_vector)