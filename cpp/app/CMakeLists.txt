project(sc-runner)

set(EXE_SOURCES
    ${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME} ${EXE_SOURCES})

target_link_libraries (${PROJECT_NAME} PRIVATE
    stratton-chu
    fftw3
)

target_include_directories(${PROJECT_NAME} PRIVATE ../stratton-chu-library/include)

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
