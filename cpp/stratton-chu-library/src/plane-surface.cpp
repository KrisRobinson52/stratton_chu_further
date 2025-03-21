#include "stratton-chu/plane-surface.hpp"


PlaneSurface::PlaneSurface(const Position& r0, const Vector& alpha, const Vector& beta) :
    m_r0(r0), m_alpha(alpha), m_beta(beta)
{
    m_dS_over_dxdy_const = m_alpha % m_beta;
}

Position PlaneSurface::point(const Vector2D& pos) const
{
    Position result;
    result = m_r0 + m_alpha * pos[0] + m_beta * pos[1];
    return result;
}

Vector2D PlaneSurface::point2d(const Vector& pos) const
{
    Vector2D result;
    Vector support=pos-m_r0;
    result[0] = support*m_alpha;
    result[1] = support*m_beta;
    return result;
}

Vector PlaneSurface::dS_over_dxdy(const Vector2D&) const
{
    return m_dS_over_dxdy_const;
}

Vector PlaneSurface::tau1(const Vector2D& pos) const
{
    return m_alpha;
}

Vector PlaneSurface::tau2(const Vector2D& pos) const
{
    return m_beta;
}
