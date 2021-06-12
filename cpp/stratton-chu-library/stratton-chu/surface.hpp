#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "stratton-chu/types.hpp"

class ISurface
{
public:
    virtual Position point(const Vector2D& pos) const = 0;
    virtual Vector dS_over_dxdy(const Vector2D& pos) const = 0;
    virtual Vector tau1(const Vector2D& pos) const = 0;
    virtual Vector tau2(const Vector2D& pos) const = 0;
    virtual ~ISurface() = default;
};

class SurfaceBase : public ISurface
{
public:
    Vector dS_over_dxdy(const Vector2D& pos) const override
    {
        return tau1(pos) % tau2(pos);
    }
};

struct SurfaceRegion
{
    SurfaceRegion(double x_min = -1.0, double x_max = 1.0, double y_min = -1.0, double y_max = 1.0) :
        x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max)
    { }

    double width() const {return x_max - x_min; }
    double height() const {return y_max - y_min; }

    double x_min, x_max;
    double y_min, y_max;
};

SurfaceRegion unite_regions(const SurfaceRegion& reg1, const SurfaceRegion& reg2);

#endif // SURFACE_HPP
