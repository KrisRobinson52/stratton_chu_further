#ifndef PARABOLICSURFACE_HPP
#define PARABOLICSURFACE_HPP

#include "stratton-chu/surface.hpp"

class ParabolicSurface : public SurfaceBase
{
public:
    ParabolicSurface (const Position& r0, const Vector& alpha, const Vector& beta, const double a, const double b);
    // Constructor of Parabolic surface r = r0 + α * x + β * y + [α×β] / |[α×β]| * (x²/a + y²/b)

    Position point(const Vector2D& pos) const override;
    //Vector dS_over_dxdy(const Vector2D& pos) const override;
    Vector tau1(const Vector2D& pos) const override;
    Vector tau2(const Vector2D& pos) const override;

private:
    Position m_r0;
    Vector m_alpha, m_beta, m_n;
    double m_a, m_b;
};

/*
class ParabolicSurfaceDistorted : public ISurface
{
public:
    ParabolicSurfaceDistorted (const Position& r0, const Vector& alpha, const Vector& beta, const double a, const double b);
    // Constructor of Parabolic surface r = r0 + α * x + β * y + [α×β] / |[α×β]| * (x²/a + y²/b)

    Position point(const Vector2D& pos) const override;
    Vector dS_over_dxdy(const Vector2D& pos) const override;
    Vector tau1(const Vector2D& pos) const override;
    Vector tau2(const Vector2D& pos) const override;

private:
    Position m_r0;
    Vector m_alpha, m_beta, m_n;
    double m_a, m_b;
};*/

#endif // PARABOLICSURFACE_HPP
