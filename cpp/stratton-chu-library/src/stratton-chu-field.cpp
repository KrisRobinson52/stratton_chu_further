#include "stratton-chu/stratton-chu-field.hpp"

#include "stratton-chu/utils.hpp"

using namespace std::complex_literals;

StrattonChuReflection::StrattonChuReflection(ISurface& surf, IField& field, const SurfaceRegion& region) :
    FieldBase(field.lambda()),
    m_surf(surf), m_field(field), m_region(region)
{

}

FieldValue StrattonChuReflection::get(const Position& pos)
{
    FieldValue result;

    result.E = integrate<VectorComplex>(
        [this, &pos] (double x, double y) { return subint_E(x, y, pos); },
        m_region.x_min, m_region.x_max,
        m_region.y_min, m_region.y_max,
        200, 200
    );
    result.E *= 1 / (4 * M_PI);
    return result;
}

VectorComplex StrattonChuReflection::subint_E(double x, double y, Vector r)
{
    Vector2D xy = {x, y};
    Vector r0 = m_surf.point(xy);
    VectorComplex N = m_surf.dS_over_dxdy(xy).cast<Complex>();

    FieldValue field = m_field.get(r0);
    GreenFunc green(r, r0, m_k);

    VectorComplex first_term = (N % field.B);
    first_term *= 2i * m_k * green.value();

    VectorComplex second_term = green.gradient();
    second_term *= 2.0 * (N * field.E);

    return first_term + second_term;
}
