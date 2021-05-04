#include "stratton-chu/utils.hpp"
#include "cubature.h"
#include "polynomial-root-finder.hpp"
#include "stratton-chu/geometry.hpp"

#include <cmath>
#include <stdexcept>

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

/*
double get_F_by_beam_parameters_p(double p, double d, double phi)
{

    double coefs[5];
    coefs[0] = sqr(p) * sqr(p+d) * tan(phi);
    coefs[1] = - 4 * p * d * (p+d);
    coefs[2] = 4 * (2*sqr(p) + 2*p*d - sqr(d)) * tan(phi);
    coefs[3] = - 16 * d;
    coefs[4] = 16 * tan(phi);

    double real[4], imag[4];
    memset(real, 0, 4 * sizeof(double));
    memset(imag, 0, 4 * sizeof(double));
    int number_of_roots = 0;

    PolynomialRootFinder finder;
    PolynomialRootFinder::RootStatus_T status = finder.FindRoots(coefs, 4, real, imag, &number_of_roots);

    //std::cout << "status = " << status;

    for (int i = 0; i < number_of_roots; i++)
    {
        if (imag[i] != 0.0)
        {
            continue;
        }

    }

    double F = real[0];
    return F;
}*/

double get_F_by_beam_parameters_alpha(double alpha, double phi, double d)
{
    double coefs[3];

    if (fabs(alpha) <= 1e-10)
    {
        coefs[0] = - sqr(d) * tan(phi);
        coefs[1] = - 4 * d;
        coefs[2] = 4 * tan(phi);
    } else {
        double ssum = tan(alpha) + tan(phi);
        double ssub = 1 - tan(alpha) * tan(phi);
        double ddiv = (1 - cos(alpha)) / sin(alpha);
        coefs[0] = - sqr(d) * ssum;
        coefs[1] = - 4 * d * (ssum * ddiv + ssub);
        coefs[2] = 4 * (ssum * (1 - sqr(ddiv)) - 2 * ssub * ddiv);
    }

    if (fabs(coefs[2]) < 1e-10)
    {
        if (fabs(coefs[1]) < 1e-10)
        {
            throw std::domain_error("Cannot find focal length with given parameters");
        }

        return -coefs[0]/coefs[1];
    }

    double D = sqr(coefs[1]) - 4 * coefs[2] * coefs[0];

    if (D < 0)
        throw std::domain_error("Cannot find focal length with given parameters");

    double F1, F2, F;
    F1 = (- coefs[1] + sqrt(D)) / (2 * coefs[2]);
    F2 = (- coefs[1] - sqrt(D)) / (2 * coefs[2]);
    F = std::max(F1, F2);
    if (F <= 0.0)
        throw std::domain_error("Cannot find focal length with given parameters");
    return F;
}

double get_p_by_beam_parameters_alpha(double alpha, double F)
{
    if (fabs(alpha) <= 1e-10)
    {
        return 0.0;
    } else {
        return 2 * F * (1 - cos(alpha)) / sin(alpha);
    }
}
