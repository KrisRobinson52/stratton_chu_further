cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project("stratton-chu-beams-cpp" CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -march=native -fext-numeric-literals")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -pg -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -g -funroll-loops -ffast-math")

add_subdirectory(stratton-chu-library)
add_subdirectory(3rdparty/googletest)
add_subdirectory(app)

enable_testing()

add_subdirectory(tests)
