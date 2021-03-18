#ifndef STRATTONCHUFIELD_HPP
#define STRATTONCHUFIELD_HPP

#include "stratton-chu/types.hpp"
#include "stratton-chu/field.hpp"
#include "stratton-chu/surface.hpp"

struct SurfaceRegion
{
    double x_min = -1.0, x_max = 1.0;
    double y_min = -1.0, y_max = 1.0;
};

class StrattonChuReflection : public FieldBase
{
public:
    StrattonChuReflection(ISurface& surf, IField& field, const SurfaceRegion& region);

    FieldValue get(const Position& pos) override;

private:

    VectorComplex subint_E(double x, double y, Vector r);

    ISurface& m_surf;
    IField& m_field;
    SurfaceRegion m_region;
};


#endif // STRATTONCHUFIELD_HPP
