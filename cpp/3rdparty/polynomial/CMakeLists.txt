project(polynomial)

set(LIB_SOURCES
    ${PROJECT_SOURCE_DIR}/polynomial-root-finder.hpp
    ${PROJECT_SOURCE_DIR}/polynomial-root-finder.cpp
)

add_library(${PROJECT_NAME} ${LIB_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
