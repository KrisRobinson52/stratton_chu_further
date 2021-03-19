#include "stratton-chu/utils.hpp"

#include "gtest/gtest.h"


TEST(IntegrateTest, Operating)
{
    auto function = [](double x, double y) { return x + 2*y; };

    double result = integrate<double>(function, -1, 1, -1, 1, 10, 10);
    ASSERT_NEAR(0.0, result, 1e-3);

    result = integrate<double>(function, 0, 1, 0, 1, 100, 100);
    ASSERT_NEAR(1.5, result, 1e-3);
}
