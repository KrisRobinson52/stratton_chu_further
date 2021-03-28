#include "stratton-chu/utils.hpp"

#include <cmath>

using namespace std::complex_literals;

GreenFunc::GreenFunc(const Position& r, const Position& r0, double k)
{
    Vector delta_r = r0 - r;
    double l = delta_r.norm();
    Complex exp_value = std::exp(Complex(0.0, 1.0) * k * l);
    m_value = exp_value / l;

    Complex grad_multiplier = 1 / sqr(l) * exp_value * (1 / l - k * Complex(0.0, 1.0));
    m_gradient = delta_r.cast<Complex>() * grad_multiplier;
}

const Complex& GreenFunc::value()
{
    return m_value;
}

const VectorComplex& GreenFunc::gradient()
{
    return m_gradient;
}
