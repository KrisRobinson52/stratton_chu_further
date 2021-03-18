cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project("stratton-chu-beams-cpp" CXX)

add_subdirectory(stratton-chu-library)
add_subdirectory(3rdparty/googletest)

enable_testing()

add_subdirectory(tests)
