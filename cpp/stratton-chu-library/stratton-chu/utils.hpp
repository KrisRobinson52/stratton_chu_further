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
ValueType integrate(std::function<ValueType(double, double)> func,
                        double x_min, double x_max, double y_min, double y_max, size_t Nx, size_t Ny)
{
    // Trapezoid integration
    double ds = (x_max - x_min) * (y_max - y_min) / (Nx * Ny);
    double step_x = (x_max - x_min) / Nx;
    double step_y = (y_max - y_min) / Ny;

    ValueType result;

    // Sum over central block
    for (size_t i = 1; i != Nx; i++)
    {
        double x = x_min + step_x * i;
        for (size_t j = 1; j != Ny; j++)
        {
            double y = y_min + step_y * j;

            ValueType df = 2 * func(x, y);

            result += df;
        }
    }

    // Sum over borders block
    for (size_t i = 0; i != Nx+1; i++)
    {
        double x = x_min + step_x * i;
        double y = y_min + step_y * i;
        result += func(x, y_min) + func(x, y_max) + func(x_min, y) + func(x_max, y);
    }

    result *= ds / 2;

    return result;
}

#endif // UTILS_HPP
