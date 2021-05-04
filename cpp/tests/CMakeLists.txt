project(unit-tests)

set(EXE_SOURCES
    ${PROJECT_SOURCE_DIR}/test-beams.cpp
    ${PROJECT_SOURCE_DIR}/test-utils.cpp
    ${PROJECT_SOURCE_DIR}/test-polynomial.cpp
)

add_executable(${PROJECT_NAME} ${EXE_SOURCES})

target_link_libraries (${PROJECT_NAME} PRIVATE
    gtest
    gtest_main
    stratton-chu
)

add_test(NAME beams-testing COMMAND ${PROJECT_NAME})
