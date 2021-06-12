#ifndef DISTORTEDSURFACE_HPP
#define DISTORTEDSURFACE_HPP

#include "stratton-chu/surface.hpp"

#include <vector>

class SurfaceDistortion : public SurfaceBase
{
public:
    struct DistortionHarmonic
    {
        double ampl, kx, ky;
    };

    SurfaceDistortion(const ISurface& pure_surface, const Vector& v, const std::vector<DistortionHarmonic>& harmonics);

    Position point(const Vector2D& pos) const override;

    Vector tau1(const Vector2D& pos) const override;

    Vector tau2(const Vector2D& pos) const override;

private:
    const ISurface& m_pure_surface;
    const std::vector<DistortionHarmonic> m_harmonics;
    const Vector m_v;
};

#endif // DISTORTEDSURFACE_HPP
