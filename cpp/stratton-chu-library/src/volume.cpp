#include "stratton-chu/volume.hpp"

VolumeRegion unite_regions(const VolumeRegion& reg1, const VolumeRegion& reg2)
{
    return VolumeRegion(
                std::min(reg1.x_min, reg2.x_min), std::max(reg1.x_max, reg2.x_max),
                std::min(reg1.y_min, reg2.y_min), std::max(reg1.y_max, reg2.y_max),
                std::min(reg1.z_min, reg2.z_min), std::max(reg1.z_max, reg2.z_max));
}
