#include "stratton-chu/plane-surface.hpp"
#include "stratton-chu/parabolic-surface.hpp"
#include "stratton-chu/parallel-beam.hpp"
#include "stratton-chu/stratton-chu-field.hpp"

#include "stratton-chu/csv-saver.hpp"

#include <iostream>

int main()
{
    std::cout << "Running stratton-chu computation" << std::endl;

    PlaneSurface surface(
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, -1.0}
    );

    double F = 1.0;
    double lambda = 0.01;
    double h = 2 * sqrt(F);

    ParabolicSurface surface2(
        {0.0, 0.0, -F},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        2.0*sqrt(F), 2.0*sqrt(F)
    );

    ParallelBeamZ beam(
        lambda,
        [](double x, double y) {
            return exp(-(sqr(x) + sqr(y)) / sqr(0.5));
        },
        [](double, double) { return 0.0; },
        0.0
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

    /*
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
    */

    StrattonChuReflection reflection2(surface2, beam2, region);

    int n_points = 101;
    double z_from = 0.0 - 0.3;
    double z_to = 0.0 + 0.3;
    double x_from = 0.0 - 0.3;
    double x_to = 0.0 + 0.3;


    CSVSaver csv_saver("output.txt");


    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            Position p(x_from + (x_to - x_from) / (n_points-1) * i, 0.0, z_from + (z_to - z_from) / (n_points-1) * j);
            FieldValue val = reflection2.get(p);
            std::cout << "Point " << p.str() << ": E = " << val.E.str() << std::endl;
            csv_saver.add_point(p, val);
        }
    }

    return 0;
}
