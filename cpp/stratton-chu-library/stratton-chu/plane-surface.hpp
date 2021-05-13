#ifndef PLANESURFACE_HPP
#define PLANESURFACE_HPP

#include "stratton-chu/surface.hpp"

class PlaneSurface : public ISurface
{
public:
    PlaneSurface(const Position& r0, const Vector& alpha, const Vector& beta);

    Position point(const Vector2D& pos) const override;
    Vector dS_over_dxdy(const Vector2D& pos) const override;

private:
    Position m_r0;
    Vector m_alpha, m_beta;
    Vector m_dS_over_dxdy_const;
};

#endif // PLANESURFACE_HPP
