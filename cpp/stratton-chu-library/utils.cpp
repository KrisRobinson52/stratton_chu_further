#include "utils.hpp"

#include <cmath>

using namespace std::complex_literals;

GreenFunc::GreenFunc(const Position& r, const Position& r0, double k)
{
    Vector delta_r = r0 - r;
    double l = delta_r.norm();
    Complex exp_value = std::exp(1i * k * l);
    value = exp_value / l;

    Complex grad_multiplier = 1 / sqr(l) * exp_value * (1 / l - k * 1i);
    gradient = delta_r.cast<Complex>() * grad_multiplier;
}
