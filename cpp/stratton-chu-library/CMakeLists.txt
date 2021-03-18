project(stratton-chu)

set(LIB_SOURCES
    ${PROJECT_SOURCE_DIR}/beams.cpp
    ${PROJECT_SOURCE_DIR}/beams.hpp
    ${PROJECT_SOURCE_DIR}/geometry.hpp
    ${PROJECT_SOURCE_DIR}/types.hpp
    ${PROJECT_SOURCE_DIR}/utils.cpp
    ${PROJECT_SOURCE_DIR}/utils.hpp
)

add_library(${PROJECT_NAME} ${LIB_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
