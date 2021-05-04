#include "polynomial-root-finder.hpp"

#include "gtest/gtest.h"

#include <iostream>
#include <cmath>

TEST(PolynomialRootFinderTest, Order2)
{
    PolynomialRootFinder finder;
    const int degree = 2;
    double coefficients[degree + 1] = {-3.0, 1.0, 2.0};
    double real[degree], imag[degree];
    memset(real, 0, degree * sizeof(double));
    memset(imag, 0, degree * sizeof(double));
    int number_of_roots = 0;
    PolynomialRootFinder::RootStatus_T status = finder.FindRoots(coefficients, degree, real, imag, &number_of_roots);
    ASSERT_EQ(status, PolynomialRootFinder::SUCCESS);
    ASSERT_EQ(number_of_roots, degree);
    double true_roots[degree] = {1, -1.5};
    for (int i = 0; i < degree; i++)
    {
        ASSERT_EQ(imag[i], 0.0);
        bool found = false;
        for (int j = 0; j < degree; j++)
        {
            if (fabs(true_roots[j] - real[i]) < 1e-5)
            {
                found = true;
                break;
            }
        }
        ASSERT_TRUE(found);
        //std::cout << "real = " << real[i] << ", imag = " << imag[i] << std::endl;
    }

}

TEST(PolynomialRootFinderTest, Order5)
{
    PolynomialRootFinder finder;
    const int degree = 5;
    double coefficients[degree + 1] = {-120.0, 94.0, 51.0, -23.0, -3.0, 1.0};
    double real[degree], imag[degree];
    memset(real, 0, degree * sizeof(double));
    memset(imag, 0, degree * sizeof(double));
    int number_of_roots = 0;
    PolynomialRootFinder::RootStatus_T status = finder.FindRoots(coefficients, degree, real, imag, &number_of_roots);
    ASSERT_EQ(status, PolynomialRootFinder::SUCCESS);
    ASSERT_EQ(number_of_roots, degree);
    double true_roots[degree] = {1.0, -2.0, 3.0, -4.0, 5.0};
    for (int i = 0; i < degree; i++)
    {
        ASSERT_EQ(imag[i], 0.0);
        bool found = false;
        for (int j = 0; j < degree; j++)
        {
            if (fabs(true_roots[j] - real[i]) < 1e-5)
            {
                found = true;
                break;
            }
        }
        ASSERT_TRUE(found);
        //std::cout << "real = " << real[i] << ", imag = " << imag[i] << std::endl;
    }

}
