cmake_minimum_required(VERSION 3.9.4 FATAL_ERROR)

project("stratton-chu" CXX)

set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_INTERPROCEDURAL_OPIMIZATION TRUE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -march=native -fext-numeric-literals -Wno-error=maybe-uninitialized")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -pg -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -flto -funroll-loops -ffast-math -march=native -fext-numeric-literals")
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
# set( CMAKE_SKIP_BUILD_RPATH true )

add_library(Zernike STATIC ${CMAKE_SOURCE_DIR}/cpp/3rdparty/zernike/src/zernike_indices.cpp
                           ${CMAKE_SOURCE_DIR}/cpp/3rdparty/zernike/src/zernike_radial_poly.cpp
                           ${CMAKE_SOURCE_DIR}/cpp/3rdparty/zernike/src/zernike_poly.cpp
                           )
include_directories(${CMAKE_SOURCE_DIR}/cpp/3rdparty/zernike/include
                    ${CMAKE_SOURCE_DIR}/cpp/3rdparty/math/include)

add_subdirectory(cpp)
