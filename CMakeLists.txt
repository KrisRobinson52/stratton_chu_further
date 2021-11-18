cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project("stratton-chu" CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -march=native -fext-numeric-literals -Wno-error=maybe-uninitialized")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -pg -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops -ffast-math")

add_subdirectory(cpp)
