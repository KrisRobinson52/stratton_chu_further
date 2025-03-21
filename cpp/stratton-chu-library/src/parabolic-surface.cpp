#include "stratton-chu/parabolic-surface.hpp"

// Constructor of Parabolic surface r = r0 + α * x + β * y + [α×β] / |[α×β]| * (x²/a + y²/b)
ParabolicSurface::ParabolicSurface(const Position& r0, const Vector& alpha, const Vector& beta,
                                   const double a, const double b) :
    m_r0(r0), m_alpha(alpha), m_beta(beta), m_a(a), m_b(b)
{
    m_n = m_alpha % m_beta;
    m_n /= m_n.norm();
}

Position ParabolicSurface::point(const Vector2D& pos) const
{
    Position result;
    result = m_n;
    result *= sqr(pos[0])/m_a + sqr(pos[1])/m_b;
    result += m_r0 + m_alpha * pos[0] + m_beta * pos[1];
    return result;
}
/*
Vector ParabolicSurface::dS_over_dxdy(const Vector2D& pos) const
{
    return tau1(pos) % tau2(pos);
}*/

Vector2D ParabolicSurface::point2d(const Vector& pos) const
{
    Vector2D result;
    Vector support=pos-m_r0;
    result[0] = support*m_alpha;
    result[1] = support*m_beta;
    return result;
}

Vector ParabolicSurface::tau1(const Vector2D& pos) const
{
    Vector r_x;

    r_x = m_n;
    r_x *= 2 * pos[0] / m_a;
    r_x += m_alpha;

    return r_x;
}

Vector ParabolicSurface::tau2(const Vector2D& pos) const
{
    Vector r_y;

    r_y = m_n;
    r_y *= 2 * pos[1] / m_b;
    r_y += m_beta;

    return r_y;
}
