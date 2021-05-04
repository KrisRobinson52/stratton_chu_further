#include "stratton-chu/utils.hpp"

#include "gtest/gtest.h"
#include <iostream>

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
    SurfaceRegion reg1(-1.0, 1.0, -1.0, 1.0);

    VectorComplex result = integrate_cubature(test_integrand, reg1, 1e-2, 1e-2);
    ASSERT_NEAR(0.0, result[0].real(), 1e-3);
}

TEST(GetFbyBeamParametersTest, OperatingByAlpha)
{
    double F = get_F_by_beam_parameters_alpha(0.0, M_PI_2, 100.0);
    ASSERT_NEAR(50.0, F, 1e-5);

    //std::cout << "[";
/*    for(double alpha = 0.0; alpha<3*M_PI_4; alpha+=0.05)
    {
        std::cout << "[" << alpha << ", ";
        double F = get_F_by_beam_parameters_alpha(alpha, M_PI_2/3, 100.0);
        std::cout << F << "]," << std::endl;

    }*/
    //std::cout << "]" << std::endl;
}
