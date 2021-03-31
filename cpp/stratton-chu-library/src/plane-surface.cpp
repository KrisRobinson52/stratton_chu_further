#include "stratton-chu/plane-surface.hpp"


PlaneSurface::PlaneSurface(const Position& r0, const Vector& alpha, const Vector& beta) :
    m_r0(r0), m_alpha(alpha), m_beta(beta)
{
    m_dS_over_dxdy_const = m_alpha % m_beta;
}

Position PlaneSurface::point(const Vector2D& pos)
{
    Position result;
    result = m_r0 + m_alpha * pos[0] + m_beta * pos[1];
    return result;
}

Vector PlaneSurface::dS_over_dxdy(const Vector2D&)
{
    return m_dS_over_dxdy_const;
}
