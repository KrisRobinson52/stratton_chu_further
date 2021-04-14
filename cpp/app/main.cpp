#include "stratton-chu/plane-surface.hpp"
#include "stratton-chu/parabolic-surface.hpp"
#include "stratton-chu/parallel-beam.hpp"
#include "stratton-chu/stratton-chu-field.hpp"

#include "stratton-chu/csv-saver.hpp"
#include "stratton-chu/vtk-saver.hpp"

#include <iostream>

int main()
{   
    /*
    VTKVolumeSaver saver(100, 100, 100, "Test quantity");

    for (size_t x = 0; x < 100; x++)
    {
        for (size_t y = 0; y < 100; y++)
        {
            for (size_t z = 0; z < 100; z++)
            {
                saver.set_point(x, y, z, Position(x, y, z), sin(y), cos(z), sin(x/10.0)*cos(y/10.0)*sin(z/10.0));
            }
        }
    }

    saver.save("test-data");*/

    /*
    VTKSurfaceSaver saver(10, 10, "Test quantity");

    for (size_t x = 0; x < 10; x++)
    {
        for (size_t y = 0; y < 10; y++)
        {
            saver.set_point(x, y, Position(x, y, 0.1*x*y), sin(y), cos(y), sin(x)*cos(y));
        }
    }

    saver.save("test-data");
*/

    std::cout << "Running stratton-chu computation" << std::endl;

    PlaneSurface surface(
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, -1.0}
    );

    ParallelBeamZ beam(
        0.1,
        [](double x, double y) {
            return exp(-(sqr(x) + sqr(y)) / sqr(0.5));
        },
        [](double, double) { return 0.0; },
        0.0
    );

    SurfaceRegion region;

    StrattonChuReflection reflection(surface, beam, region);

    int n_points = 10;

    double z_from = -2.0;
    double z_to = 2.0;

    CSVSaver csv_saver("output.txt");

    for (int i = 0; i < n_points; i++)
    {
        Position p(0.0, 2.05, z_from + (z_to - z_from) / (n_points-1) * i);
        FieldValue val = reflection.get(p);
        std::cout << "Point " << p.str() << ": E = " << val.E.str() << std::endl;
        csv_saver.add_point(p, val);
    }

    return 0;
}
