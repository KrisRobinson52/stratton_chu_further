project(stratton-chu)

find_package(VTK COMPONENTS
  vtkCommonColor
  vtkCommonCore
  vtkCommonDataModel
  vtkFiltersGeneral
  vtkIOXML
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  REQUIRED
)

set(LIB_SOURCES
    ${PROJECT_SOURCE_DIR}/stratton-chu/parallel-beam.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/plane-surface.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/stratton-chu-field.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/field.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/surface.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/geometry.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/types.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/utils.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/parabolic-surface.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/csv-saver.hpp
    ${PROJECT_SOURCE_DIR}/stratton-chu/vtk-saver.hpp

    ${PROJECT_SOURCE_DIR}/src/parallel-beam.cpp
    ${PROJECT_SOURCE_DIR}/src/utils.cpp
    ${PROJECT_SOURCE_DIR}/src/field.cpp
    ${PROJECT_SOURCE_DIR}/src/stratton-chu-field.cpp
    ${PROJECT_SOURCE_DIR}/src/plane-surface.cpp
    ${PROJECT_SOURCE_DIR}/src/parabolic-surface.cpp
    ${PROJECT_SOURCE_DIR}/src/csv-saver.cpp
    ${PROJECT_SOURCE_DIR}/src/vtk-saver.cpp
)

add_library(${PROJECT_NAME} ${LIB_SOURCES})

include(${VTK_USE_FILE})

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC ${VTK_LIBRARIES} cubature)

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
