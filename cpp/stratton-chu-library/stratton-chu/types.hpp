#ifndef TYPES_HPP
#define TYPES_HPP

#include "stratton-chu/geometry.hpp"

#include <complex>

using Complex = std::complex<double>;
using Position = StaticVector<3>;
using Vector = StaticVector<3>;
using Vector2D = StaticVector<2>;
using VectorComplex = StaticVector<3, Complex>;

struct SurfaceRegion
{
    SurfaceRegion(double x_min = -1.0, double x_max = 1.0, double y_min = -1.0, double y_max = 1.0) :
        x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max)
    { }

    double x_min = -1.0, x_max = 1.0;
    double y_min = -1.0, y_max = 1.0;
};

constexpr double c = 29979245800.0;

#endif // TYPES_HPP
