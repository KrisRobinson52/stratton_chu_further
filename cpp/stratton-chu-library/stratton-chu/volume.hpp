#ifndef VOLUME_HPP
#define VOLUME_HPP

#include "stratton-chu/types.hpp"

class IVolume
{
public:
    virtual Position point(const Vector& pos) const = 0;
    virtual Vector tau1(const Vector& pos) const = 0;
    virtual Vector tau2(const Vector& pos) const = 0;
    virtual Vector tau3(const Vector& pos) const = 0;
    virtual ~IVolume() = default;
};

struct VolumeRegion
{
    VolumeRegion(double x_min = -1.0, double x_max = 1.0, double y_min = -1.0, double y_max = 1.0, double z_min = -1.0, double z_max = 1.0) :
        x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), z_min(z_min), z_max(z_max)
    { }

    double widthx() const {return x_max - x_min; }
    double widthy() const {return y_max - y_min; }
    double height() const {return z_max - z_min; }

    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
};

VolumeRegion unite_regions(const VolumeRegion& reg1, const VolumeRegion& reg2);

#endif // VOLUME_HPP
