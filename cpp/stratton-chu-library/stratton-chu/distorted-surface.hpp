#ifndef DISTORTEDSURFACE_HPP
#define DISTORTEDSURFACE_HPP

#include "stratton-chu/surface.hpp"

#include <vector>

class SurfaceDistortionHarmonic : public SurfaceBase
{
public:
    struct DistortionHarmonic
    {
        double ampl, kx, ky;
    };

    SurfaceDistortionHarmonic(const ISurface& pure_surface, const Vector& v, const std::vector<DistortionHarmonic>& harmonics);

    Position point(const Vector2D& pos) const override;

    Vector tau1(const Vector2D& pos) const override;

    Vector tau2(const Vector2D& pos) const override;

private:
    const ISurface& m_pure_surface;
    const std::vector<DistortionHarmonic> m_harmonics;
    const Vector m_v;
};

class SurfaceDistortionLegendre : public SurfaceBase
{
public:
    struct DistortionPolinom
    {
        DistortionPolinom(double ampl, double alpha, int number);

        const double ampl;
        const Vector2D direction;
        const int number;
    };

    SurfaceDistortionLegendre(const ISurface& pure_surface, const Vector& v, double radius, Vector2D center, const std::vector<DistortionPolinom>& harmonics);

    Position point(const Vector2D& pos) const override;

    Vector tau1(const Vector2D& pos) const override;

    Vector tau2(const Vector2D& pos) const override;

private:
    const ISurface& m_pure_surface;
    const std::vector<DistortionPolinom> m_harmonics;
    const Vector m_v;
    const Vector2D m_center;
    const double m_radius;
};


#endif // DISTORTEDSURFACE_HPP
