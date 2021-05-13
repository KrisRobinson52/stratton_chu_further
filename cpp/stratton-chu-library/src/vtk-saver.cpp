#include "stratton-chu/vtk-saver.hpp"

#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>

VTKSurfaceSaver::VTKSurfaceSaver(size_t Nx, size_t Ny, const char* name) :
    m_Nx(Nx), m_Ny(Ny)
{
    size_t points_count = m_Nx * m_Ny;
    m_points->SetNumberOfPoints(points_count);

    m_data->SetName(name);
    m_data->SetNumberOfComponents(3);
    m_data->SetComponentName(0, "X");
    m_data->SetComponentName(1, "Y");
    m_data->SetComponentName(2, "Z");
    m_data->SetNumberOfTuples(points_count);

    build_surface();
}

void VTKSurfaceSaver::set_point(size_t x_index, size_t y_index, Position p, double value_x, double value_y, double value_z)
{
    size_t glob_ind = global_index(x_index, y_index);
    m_points->SetPoint(glob_ind, p.x);

    double point_data[] = {value_x, value_y, value_z};
    m_data->SetTuple(glob_ind, point_data);
}

void VTKSurfaceSaver::set_point(size_t x_index, size_t y_index, Position p, Vector vec)
{
    set_point(x_index, y_index, p, vec[0], vec[1], vec[2]);
}

void VTKSurfaceSaver::save(const char* filename)
{
    std::string full_name(filename);
    full_name += ".vtp";
    vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetFileName(full_name.c_str());
    writer->SetInputData(m_polydata);
    writer->Write();
}

size_t VTKSurfaceSaver::global_index(size_t x, size_t y)
{
    return x * m_Ny + y;
}

void VTKSurfaceSaver::build_surface()
{
    for (size_t x = 0; x < m_Nx-1; x++)
    {
        for (size_t y = 0; y < m_Ny-1; y++)
        {
            vtkNew<vtkTriangle> t1;
            t1->GetPointIds()->SetId(0, global_index(x, y));
            t1->GetPointIds()->SetId(1, global_index(x+1, y));
            t1->GetPointIds()->SetId(2, global_index(x, y+1));

            vtkNew<vtkTriangle> t2;
            t2->GetPointIds()->SetId(0, global_index(x+1, y));
            t2->GetPointIds()->SetId(1, global_index(x+1, y+1));
            t2->GetPointIds()->SetId(2, global_index(x, y+1));

            m_triangles->InsertNextCell(t1);
            m_triangles->InsertNextCell(t2);
        }
    }

    m_polydata->SetPoints(m_points);
    m_polydata->SetPolys(m_triangles);
    m_polydata->GetPointData()->SetScalars(m_data);
}

VTKVolumeSaver::VTKVolumeSaver(size_t Nx, size_t Ny, size_t Nz, const char* name) :
    m_Nx(Nx), m_Ny(Ny), m_Nz(Nz)
{
    size_t points_count = m_Nx * m_Ny * m_Nz;
    m_points->SetNumberOfPoints(points_count);

    m_data->SetName(name);
    m_data->SetNumberOfComponents(3);
    m_data->SetComponentName(0, "X");
    m_data->SetComponentName(1, "Y");
    m_data->SetComponentName(2, "Z");
    m_data->SetNumberOfTuples(points_count);

    build_volume();
}

void VTKVolumeSaver::set_point(size_t x_index, size_t y_index, size_t z_index, Position p, double value_x, double value_y, double value_z)
{
    size_t glob_ind = global_index(x_index, y_index, z_index);
    m_points->SetPoint(glob_ind, p.x);

    double point_data[] = {value_x, value_y, value_z};
    m_data->SetTuple(glob_ind, point_data);
}

size_t VTKVolumeSaver::global_index(size_t x, size_t y, size_t z)
{
    return x * m_Ny * m_Nz + y * m_Nz + z;
}

void VTKVolumeSaver::build_volume()
{
    for (size_t x = 0; x < m_Nx-1; x++)
    {
        for (size_t y = 0; y < m_Ny-1; y++)
        {
            for (size_t z = 0; z < m_Nz-1; z++)
            {
                vtkNew<vtkHexahedron> hex;
                hex->GetPointIds()->SetId(0, global_index(x,     y, z));
                hex->GetPointIds()->SetId(1, global_index(x+1,   y, z));
                hex->GetPointIds()->SetId(2, global_index(x+1, y+1, z));
                hex->GetPointIds()->SetId(3, global_index(x,   y+1, z));

                hex->GetPointIds()->SetId(4, global_index(x,     y, z+1));
                hex->GetPointIds()->SetId(5, global_index(x+1,   y, z+1));
                hex->GetPointIds()->SetId(6, global_index(x+1, y+1, z+1));
                hex->GetPointIds()->SetId(7, global_index(x,   y+1, z+1));


                m_hexagons->InsertNextCell(hex);
                m_grid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
            }
        }
    }

    m_grid->SetPoints(m_points);
    m_grid->GetPointData()->SetScalars(m_data);
}

void VTKVolumeSaver::save(const char* filename)
{
    std::string full_name(filename);
    full_name += ".vtu";
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(full_name.c_str());
    writer->SetInputData(m_grid);
    writer->Write();
}
