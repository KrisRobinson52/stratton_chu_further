#ifndef PARABOLICSURFACE_HPP
#define PARABOLICSURFACE_HPP

#include "stratton-chu/surface.hpp"

class ParabolicSurface : public ISurface
{
public:
    ParabolicSurface (const Position& r0, const Vector& alpha, const Vector& beta, const double a, const double b);
    // Constructor of Parabolic surface r = r0 + α * x + β * y + [α×β] / |[α×β]| * (x²/a² + y²/b²)

    Position point(const Vector2D& pos) override;
    Vector dS_over_dxdy(const Vector2D& pos) override;

private:
    Position m_r0;
    Vector m_alpha, m_beta, m_n;
    double m_a, m_b;
};

#endif // PARABOLICSURFACE_HPP
