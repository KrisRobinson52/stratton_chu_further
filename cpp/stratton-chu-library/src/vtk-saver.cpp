#include "stratton-chu/vtk-saver.hpp"

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkXMLPolyDataReader.h>
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

void VTKSurfaceSaver::save(const char* filename)
{
    vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetFileName(filename);
    writer->SetInputData(m_polydata);
    writer->Write();
}

size_t VTKSurfaceSaver::global_index(size_t x, size_t y)
{
    return x + y * m_Nx;
}

void VTKSurfaceSaver::build_surface()
{
    for (size_t x = 0; x < m_Nx-1; x++)
    {
        for (size_t y = 0; y < m_Nx-1; y++)
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
