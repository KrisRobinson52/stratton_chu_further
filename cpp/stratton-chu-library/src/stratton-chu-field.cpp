#include "stratton-chu/stratton-chu-field.hpp"

#include "stratton-chu/utils.hpp"

using namespace std::complex_literals;

StrattonChuReflection::StrattonChuReflection(ISurface& surf, IField& field, const SurfaceRegion& region) :
    FieldBase(field.lambda()),
    m_surf(surf), m_field(field), m_region(region)
{
}

FieldValue StrattonChuReflection::get(const Position& pos) const
{
    FieldValue result;
/*
    result.E = integrate_trapezoid<VectorComplex>(
        [this, &pos] (double x, double y) -> VectorComplex { return subint_E(x, y, pos); },
        m_region,
        200, 200
    );*/

    result.E = integrate_cubature(
        [this, &pos] (double x, double y) -> VectorComplex { return subint_E(x, y, pos); },
        m_region,
        1e-4, 1e-4
    );
    result.E *= 1 / (4 * M_PI);

    // @TODO Calculate B

    return result;
}

VectorComplex StrattonChuReflection::subint_E(double x, double y, const Position& r) const
{
    Vector2D xy = {x, y};
    Vector r0 = m_surf.point(xy);
    VectorComplex N = m_surf.dS_over_dxdy(xy).cast<Complex>();

    FieldValue field = m_field.get(r0);
    GreenFunc green(r, r0, m_k);

    VectorComplex first_term = (N % field.B);
    first_term *= 2.0 * Complex(0.0, 1.0) * m_k * green.value();

    VectorComplex second_term = green.gradient();
    second_term *= 2.0 * (N * field.E);

    return first_term + second_term;
}

VectorComplex StrattonChuReflection::subint_B(double x, double y, const Position& r) const
{
    /*Vector2D xy = {x, y};
    Vector r0 = m_surf.point(xy);
    VectorComplex N = m_surf.dS_over_dxdy(xy).cast<Complex>();

    FieldValue field = m_field.get(r0);
    GreenFunc green(r, r0, m_k);

    VectorComplex first_term = (N % field.B);
    first_term *= 2.0 * Complex(0.0, 1.0) * m_k * green.value();

    VectorComplex second_term = green.gradient();
    second_term *= 2.0 * (N * field.E);

    return first_term + second_term;
    */
    return VectorComplex();
}
