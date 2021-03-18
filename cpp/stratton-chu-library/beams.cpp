#include "beams.hpp"
#include "utils.hpp"

#include <cmath>

using namespace std::complex_literals;

PlaneSurface::PlaneSurface(const StaticVector<3>& r0, const StaticVector<3>& alpha, const StaticVector<3>& beta) :
    r0(r0), alpha(alpha), beta(beta)
{
    dS_over_dxdy_const = alpha % beta;
}

StaticVector<3> PlaneSurface::point(const StaticVector<2>& pos)
{
    StaticVector<3> result;
    result = r0 + alpha * pos[0] + beta * pos[1];
    return result;
}

StaticVector<3> PlaneSurface::dS_over_dxdy(const StaticVector<2>&)
{
    return dS_over_dxdy_const;
}

FieldBase::FieldBase(double lambda) :
    m_lambda(lambda)
{
    m_k = 2 * M_PI / lambda;
    m_omega = m_k * c;
}

double FieldBase::k()
{
    return m_k;
}

double FieldBase::omega()
{
    return m_omega;
}

double FieldBase::lambda()
{
    return m_lambda;
}

ParallelBeamZ::ParallelBeamZ(double lambda, FieldAmpl Ex, FieldAmpl Ey, double z0) :
    FieldBase(lambda), m_Ex(Ex), m_Ey(Ey), m_z0(z0)
{
}

FieldValue ParallelBeamZ::get(const Position& pos)
{
    double x = pos[0], y = pos[1], z = pos[2];
    std::complex phase_mul = std::exp( -1i * m_k *(z - m_z0));
    FieldValue result;

    std::complex Ex = m_Ex(x, y) * phase_mul;
    std::complex Ey = m_Ey(x, y) * phase_mul;

    result.E[0] = Ex; result.E[0] = Ey; result.E[0] = 0;
    result.B[0] = -Ey; result.B[0] = Ex; result.E[0] = 0;

    return result;
}

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
    GreenFunc green(r0, r, m_k);

    VectorComplex first_term = (N % field.B);
    first_term *= 2i * m_k * green.value;

    VectorComplex second_term = green.gradient;
    second_term *= 2.0 * (N * field.E);

    return first_term + second_term;
}
