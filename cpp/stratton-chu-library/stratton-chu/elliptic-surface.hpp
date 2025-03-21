#ifndef ELLIPTICSURFACE_HPP
#define ELLIPTICSURFACE_HPP

#include "stratton-chu/surface.hpp"

class EllipticSurface : public SurfaceBase
{
public:
    EllipticSurface (const Position& r0, const Vector& alpha, const Vector& beta, const double a, const double b, const double c);
    // Constructor of half of Elliptic surface r = r0 + α * x + β * y + [α×β] / |[α×β]| * c*(1 - x²/a - y²/b)^0.5

    bool BelongsToSurface (const Vector2D& pos);
    
    Position point(const Vector2D& pos) const override;
    Vector2D point2d(const Vector& pos) const override;
    Vector dS_over_dxdy(const Vector2D& pos) const override;
    Vector tau1(const Vector2D& pos) const override;
    Vector tau2(const Vector2D& pos) const override;

private:
    Position m_r0;
    Vector m_alpha, m_beta, m_n;
    double m_a, m_b, m_c;
};


#endif // ELLIPTICSURFACE_HPP