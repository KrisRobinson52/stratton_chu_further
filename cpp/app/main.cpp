#include "stratton-chu/plane-surface.hpp"
#include "stratton-chu/parabolic-surface.hpp"
#include "stratton-chu/parallel-beam.hpp"
#include "stratton-chu/stratton-chu-field.hpp"

#include "stratton-chu/csv-saver.hpp"
#include "stratton-chu/vtk-saver.hpp"

#include <tbb/tbb.h>
#include <iostream>

int main()
{   

    std::cout << "Running stratton-chu computation" << std::endl;

    double F = 1.0;
    double lambda = 0.000091; // 910 nm
    double h = 2 * F;

    ParabolicSurface surface2(
        {0.0, 0.0, -F},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        4.0*F, 4.0*F
    );

    ParallelBeamZ beam2(
        lambda,
        [h](double x, double y) { return exp( - (sqr(x - h) + sqr(y)) / sqr(0.5) / 2); },
        [](double, double) { return 0.0; },
        0.0
    );

    SurfaceRegion region;
    region.x_min = 0.0;
    region.x_max = 4.0;
    region.y_min = -2.0;
    region.y_max = 2.0;

    StrattonChuReflection reflection2(surface2, beam2, region);

    int n_points = 201;
    double z_from = 0.0 - 0.001;
    double z_to = 0.0 + 0.001;
    double x_from = 0.0 - 0.001;
    double x_to = 0.0 + 0.001;

    //CSVSaver csv_saver("output.txt");
    VTKSurfaceSaver vtk_saver(n_points, n_points, "Field E");

    tbb::parallel_for( int(0), n_points,
        [x_from, x_to, z_from, z_to, n_points, &reflection2, &vtk_saver](int i) {
            for (int j = 0; j < n_points; j++)
            {
                Position p(x_from + (x_to - x_from) / (n_points-1) * i, 0.0, z_from + (z_to - z_from) / (n_points-1) * j);
                FieldValue val = reflection2.get(p);
                std::cout << "Point " << p.str() << ": E = " << val.E.str() << std::endl;
                //csv_saver.add_point(p, val);
                vtk_saver.set_point(i, j, p, val.E.x[0].real(), val.E.x[1].real(), val.E.x[2].real());
            }
        }
    );

    /*
    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            Position p(x_from + (x_to - x_from) / (n_points-1) * i, 0.0, z_from + (z_to - z_from) / (n_points-1) * j);
            FieldValue val = reflection2.get(p);
            std::cout << "Point " << p.str() << ": E = " << val.E.str() << std::endl;
            csv_saver.add_point(p, val);
            vtk_saver.set_point(i, j, p, val.E.x[0].real(), val.E.x[1].real(), val.E.x[2].real());
        }
    }*/
    vtk_saver.save("NU_SHTO");

    return 0;
}
