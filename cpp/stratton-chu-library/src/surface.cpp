#include "stratton-chu/surface.hpp"

SurfaceRegion unite_regions(const SurfaceRegion& reg1, const SurfaceRegion& reg2)
{
    return SurfaceRegion(
                std::min(reg1.x_min, reg2.x_min), std::max(reg1.x_max, reg2.x_max),
                std::min(reg1.y_min, reg2.y_min), std::max(reg1.y_max, reg2.y_max));
}
