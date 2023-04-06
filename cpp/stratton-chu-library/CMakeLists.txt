project(stratton-chu)

# find_package(VTK COMPONENTS
#   CommonCore
#   CommonDataModel
#   IOXML
#   REQUIRED
# )

set(LIB_SOURCES
    ${PROJECT_SOURCE_DIR}/stratton-chu/utils.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/types.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/geometry.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/surface.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/plane-surface.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/parabolic-surface.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/parallelepiped-volume.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/distorted-surface.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/field.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/parallel-beam.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/stratton-chu-field.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/csv-saver.hpp
    # ${PROJECT_SOURCE_DIR}/stratton-chu/vtk-saver.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/spec-inverse-fourier.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/fields-in-region.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/volume.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/fields-in-3dregion.hpp


    ${PROJECT_SOURCE_DIR}/src/utils.cpp
    ${PROJECT_SOURCE_DIR}/src/surface.cpp
    ${PROJECT_SOURCE_DIR}/src/plane-surface.cpp
    ${PROJECT_SOURCE_DIR}/src/parabolic-surface.cpp
    ${PROJECT_SOURCE_DIR}/src/parallelepiped-volume.cpp
    ${PROJECT_SOURCE_DIR}/src/distorted-surface.cpp
    ${PROJECT_SOURCE_DIR}/src/field.cpp
    ${PROJECT_SOURCE_DIR}/src/parallel-beam.cpp
    ${PROJECT_SOURCE_DIR}/src/stratton-chu-field.cpp
    ${PROJECT_SOURCE_DIR}/src/csv-saver.cpp
    # ${PROJECT_SOURCE_DIR}/src/vtk-saver.cpp
    ${PROJECT_SOURCE_DIR}/src/spec-inverse-fourier.cpp
    ${PROJECT_SOURCE_DIR}/src/volume.cpp

)

add_library(${PROJECT_NAME} ${LIB_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../3rdparty/cache)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../3rdparty/tracker)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../3rdparty/threads)

target_link_libraries(${PROJECT_NAME} PUBLIC
  # VTK::CommonCore
  # VTK::CommonDataModel
  # VTK::IOXML
  cubature
  polynomial
  boost_math
)

