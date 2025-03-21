project(sc-runner)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
#find_package(ALGLIB REQUIRED)

set(EXE_SOURCES
    ${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME} ${EXE_SOURCES})

target_link_libraries (${PROJECT_NAME} PRIVATE
    stratton-chu
    # alglib
    #fftw3
    Threads::Threads
)

target_include_directories(${PROJECT_NAME} PRIVATE ../stratton-chu-library/include)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../3rdparty/tracker)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../3rdparty/cache)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../3rdparty/threads)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../3rdparty/zernike)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../3rdparty/alglib/src/alglib/src)