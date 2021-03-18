#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "stratton-chu/types.hpp"

class ISurface
{
public:
    virtual Position point(const Vector2D& pos) = 0;
    virtual Vector dS_over_dxdy(const Vector2D& pos) = 0;
    virtual ~ISurface() = default;
};

#endif // SURFACE_HPP
