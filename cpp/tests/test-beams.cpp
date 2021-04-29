#include "stratton-chu/plane-surface.hpp"
#include "stratton-chu/parallel-beam.hpp"
#include "stratton-chu/stratton-chu-field.hpp"

#include "gtest/gtest.h"

#include <cmath>

TEST(PlaneSurfaceTest, Operating)
{
    PlaneSurface surf({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    auto p = surf.point({1, 2});
    ASSERT_NEAR(p[0], 1.0, 1e-5);
    ASSERT_NEAR(p[1], 2.0, 1e-5);
    ASSERT_NEAR(p[2], 0.0, 1e-5);
}

TEST(StrattonChuTest, PlaneSurfReflection)
{
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

    for (int i = 0; i < n_points; i++)
    {
        Position p(0.0, 2.05, z_from + (z_to - z_from) / (n_points-1) * i);
        FieldValue val = reflection.get(p);
        std::cout << "Point " << p.str() << ": " << val.E.str() << std::endl;
    }
}

TEST(ParallelBeamAlphaTest, CompareWithParallelBeamZ)
{
    FieldAmpl E1_ampl_z = [](double x, double y) {
        return exp(-(sqr(x) + sqr(y)) / sqr(0.5));
    };

    FieldAmpl E1_ampl_a = [E1_ampl_z](double x, double y) {
        return E1_ampl_z(y, x);
    };


    FieldAmpl E2_ampl_z = [](double x, double y) {
        return exp(-(sqr(2*x) + sqr(3*y)) / sqr(1.5));
    };

    FieldAmpl E2_ampl_a = [E2_ampl_z](double x, double y) {
        return E2_ampl_z(y, x);
    };


    ParallelBeamZ beam_z(0.1, E1_ampl_z, E2_ampl_z, 0.0);

    ParallelBeamAlpha beam_a(0.1, Position(0.0, 0.0, 0.0), Vector(0.0, 3.0, 0.0), Vector(1.0, 0.0, 0.0),
                             E2_ampl_a, E1_ampl_a);

    std::vector<Position> points;

    points.push_back(Position (0.0, 0.0, 0.0));
    points.push_back(Position (3.3, 2.2, -5.5));
    points.push_back(Position (-1.0, 0.0, 0.0));

    for (size_t i = 0; i < points.size(); i++) {
        FieldValue b_z = beam_z.get(points[i]);
        FieldValue b_a = beam_a.get(points[i]);
        VectorComplex delta_E = b_z.E - b_a.E;
        VectorComplex delta_B = b_z.B - b_a.B;
        ASSERT_LE( delta_E.norm() , 1e-5) << "For point " << i;
        ASSERT_LE( delta_B.norm() , 1e-5) << "For point " << i;
    };

}















