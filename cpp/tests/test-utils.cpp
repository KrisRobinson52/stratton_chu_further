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

TEST(MaxFieldTest, Linear)
{
    VectorComplex E(1.0, 2.0, 3.0);
    E *= std::exp( Complex(0.0, 1.0) * 0.6);
    Vector res = max_field(E);
    ASSERT_NEAR(res[0], 1.0, 1e-4);
    ASSERT_NEAR(res[1], 2.0, 1e-4);
    ASSERT_NEAR(res[2], 3.0, 1e-4);
}

TEST(MaxFieldTest, Elliptic)
{
    VectorComplex E(1.0, 2.0, 0.0);
    E[1] *= std::exp( Complex(0.0, 1.0) * M_PI_2);
    Vector res = max_field(E);
    ASSERT_NEAR(res[0], 0.0, 1e-4);
    ASSERT_NEAR(res[1], 2.0, 1e-4);
    ASSERT_NEAR(res[2], 0.0, 1e-4);
}

TEST(MaxFieldTest, Circular)
{
    VectorComplex E(2.0, 2.0, 0.0);
    E[1] *= std::exp( Complex(0.0, 1.0) * M_PI_2);
    Vector res = max_field(E);
    ASSERT_NEAR(res.norm(), 2.0, 1e-4);
}

TEST(GreenFunc, Differentiating)
{
    double epsilon = 1e-8;
    Position r(-100.0, 100.0, 200.0);
    double k = 1000;
    GreenFunc gf(r, Position(1.0, 4.0, 9.0), k);
    auto val = gf.value();
    auto grad = gf.gradient();

    GreenFunc gfx(r, Position(1.0+epsilon, 4.0, 9.0), k);
    GreenFunc gfy(r, Position(1.0, 4.0+epsilon, 9.0), k);
    GreenFunc gfz(r, Position(1.0, 4.0, 9.0+epsilon), k);

    double delta_x = sqrt(std::norm(grad[0] - (gfx.value() - val) / epsilon));
    double delta_y = sqrt(std::norm(grad[1] - (gfy.value() - val) / epsilon));
    double delta_z = sqrt(std::norm(grad[2] - (gfz.value() - val) / epsilon));

    ASSERT_LE(delta_x, 1e-4);
    ASSERT_LE(delta_y, 1e-4);
    ASSERT_LE(delta_z, 1e-4);

}
