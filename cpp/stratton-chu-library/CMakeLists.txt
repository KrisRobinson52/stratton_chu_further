project(stratton-chu)

set(LIB_SOURCES
    ${PROJECT_SOURCE_DIR}/stratton-chu/parallel-beam.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/plane-surface.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/stratton-chu-field.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/field.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/surface.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/geometry.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/types.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/utils.hpp

    ${PROJECT_SOURCE_DIR}/src/parallel-beam.cpp
    ${PROJECT_SOURCE_DIR}/src/utils.cpp
    ${PROJECT_SOURCE_DIR}/src/field.cpp
    ${PROJECT_SOURCE_DIR}/src/stratton-chu-field.cpp
    ${PROJECT_SOURCE_DIR}/src/plane-surface.cpp
)

add_library(${PROJECT_NAME} ${LIB_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
