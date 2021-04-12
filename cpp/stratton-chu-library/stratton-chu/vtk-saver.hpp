#ifndef VTKSAVER_HPP
#define VTKSAVER_HPP

#include "stratton-chu/types.hpp"
#include "stratton-chu/field.hpp"

#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>

#include <vtkNew.h>

#include <string>

class VTKSurfaceSaver
{
public:
    VTKSurfaceSaver(size_t Nx, size_t Ny, const char* name);

    void set_point(size_t x_index, size_t y_index, Position p, double value_x, double value_y, double value_z);
    void save(const char* filename);

private:
    size_t global_index(size_t x, size_t y);
    void build_surface();

    size_t m_Nx, m_Ny;

    vtkNew<vtkPoints> m_points;
    vtkNew<vtkDoubleArray> m_data;
    vtkNew<vtkCellArray> m_triangles;
    vtkNew<vtkPolyData> m_polydata;
};



#endif // VTKSAVER_HPP
