#include "stratton-chu/utils.hpp"

#include "gtest/gtest.h"


TEST(IntegrateTest, Operating)
{
    auto function = [](double x, double y) { return x + 2*y; };

    SurfaceRegion reg1(-1.0, 1.0, -1.0, 1.0);
    double result = integrate_trapezoid<double>(function, reg1, 10, 10);
    ASSERT_NEAR(0.0, result, 1e-3);

    SurfaceRegion reg2(0, 1.0, 0, 1.0);
    result = integrate_trapezoid<double>(function, reg2, 100, 100);
    ASSERT_NEAR(1.5, result, 1e-3);
}

VectorComplex test_integrand(double x, double y)
{
    VectorComplex result;
    result[0].real(x + 2*y);
    return result;
}

TEST(IntegrateCubature, Operating)
{
    auto function = [](double x, double y) { return x + 2*y; };

    SurfaceRegion reg1(-1.0, 1.0, -1.0, 1.0);

    VectorComplex result = integrate_cubature(test_integrand, reg1, 1e-2, 1e-2);
    ASSERT_NEAR(0.0, result[0].real(), 1e-3);
}
