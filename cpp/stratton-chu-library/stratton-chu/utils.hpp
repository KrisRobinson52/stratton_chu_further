#ifndef UTILS_HPP
#define UTILS_HPP

#include "stratton-chu/types.hpp"

#include <complex>
#include <functional>

class GreenFunc
{
public:
    GreenFunc(const Position& r, const Position& r0, double k);

    const Complex& value();
    const VectorComplex& gradient();

private:
    Complex m_value;
    VectorComplex m_gradient;
};

template<typename ValueType>
ValueType integrate_trapezoid(std::function<ValueType(double, double)> func,
                        const SurfaceRegion& region, size_t Nx, size_t Ny)
{
    // Trapezoid integration
    double ds = (region.x_max - region.x_min) * (region.y_max - region.y_min) / (Nx * Ny);
    double step_x = (region.x_max - region.x_min) / Nx;
    double step_y = (region.y_max - region.y_min) / Ny;

    ValueType result = ValueType();

    // Sum over central block
    for (size_t i = 1; i != Nx; i++)
    {
        double x = region.x_min + step_x * i;
        for (size_t j = 1; j != Ny; j++)
        {
            double y = region.y_min + step_y * j;

            ValueType df = 2 * func(x, y);

            result += df;
        }
    }

    // Sum over borders block
    for (size_t i = 0; i != Nx+1; i++)
    {
        double x = region.x_min + step_x * i;
        double y = region.y_min + step_y * i;
        result += func(x, region.y_min) + func(x, region.y_max) + func(region.x_min, y) + func(region.x_max, y);
    }

    result *= ds / 2;

    return result;
}


VectorComplex integrate_cubature(std::function<VectorComplex(double, double)> func, const SurfaceRegion& region, double rel_tol, double abs_tol);

double get_F_by_beam_parameters_p(double p, double d, double phi);

double get_F_by_beam_parameters_alpha(double alpha, double phi, double d);

double get_p_by_beam_parameters_alpha(double alpha, double F);

#endif // UTILS_HPP
