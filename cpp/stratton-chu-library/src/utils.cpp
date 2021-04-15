#include "stratton-chu/utils.hpp"
#include "cubature.h"

#include <cmath>

using namespace std::complex_literals;

extern "C" {
    static int c_integrand(unsigned ndim, const double *x, void *body, unsigned fdim, double *fval);
}

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

VectorComplex integrate_cubature(std::function<VectorComplex(double, double)> func, const SurfaceRegion& region, double rel_tol, double abs_tol)
{
    double xmin[2] = { region.x_min, region.y_min };
    double xmax[2] = { region.x_max, region.y_max };

    double result[6];
    double error[6];

    memset(result, 0, sizeof(result));
    memset(error, 0, sizeof(result));

    int ret = pcubature(6, c_integrand, reinterpret_cast<void*>(&func), 2, xmin, xmax, 0, abs_tol, rel_tol, ERROR_L2, result, error);

    VectorComplex result_vector_complex;
    result_vector_complex[0].real(result[0]);
    result_vector_complex[0].imag(result[1]);

    result_vector_complex[1].real(result[2]);
    result_vector_complex[1].imag(result[3]);

    result_vector_complex[2].real(result[4]);
    result_vector_complex[2].imag(result[5]);
    return result_vector_complex;
}

extern "C" int c_integrand(unsigned ndim, const double *x, void *body, unsigned fdim, double *fval)
{
    std::function<VectorComplex(double, double)>* func = reinterpret_cast<std::function<VectorComplex(double, double)>*>(body);
    double pos_x = x[0], pos_y = x[1];
    VectorComplex result = (*func)(pos_x, pos_y);
    fval[0] = result[0].real();
    fval[1] = result[0].imag();
    fval[2] = result[1].real();
    fval[3] = result[1].imag();
    fval[4] = result[2].real();
    fval[5] = result[2].imag();
    return 0;
}
