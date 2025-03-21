project(unit-tests)

set(EXE_SOURCES
    ${PROJECT_SOURCE_DIR}/test-beams.cpp
    ${PROJECT_SOURCE_DIR}/test-utils.cpp
    ${PROJECT_SOURCE_DIR}/test-polynomial.cpp
)

add_executable(${PROJECT_NAME} ${EXE_SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../3rdparty/alglib/src/alglib/src)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../3rdparty/threads)

target_link_libraries (${PROJECT_NAME} PRIVATE
    gtest
    gtest_main
    stratton-chu
    # alglib
)

add_test(NAME beams-testing COMMAND ${PROJECT_NAME})
