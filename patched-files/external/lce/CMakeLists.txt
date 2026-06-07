cmake_minimum_required(VERSION 3.14)
project(lce)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mcx16")

find_package(OpenMP REQUIRED)

# Automically add used submodules for the project
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
    execute_process(COMMAND "cp -rf patched-files/* ." WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

option(LCE_USE_SDSL "build lce data structure that depends on SDSL" OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/extlib")
add_subdirectory(include)
add_subdirectory(extlib)
