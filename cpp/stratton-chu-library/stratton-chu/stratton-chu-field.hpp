#ifndef STRATTONCHUFIELD_HPP
#define STRATTONCHUFIELD_HPP

#include "stratton-chu/types.hpp"
#include "stratton-chu/field.hpp"
#include "stratton-chu/surface.hpp"

class StrattonChuReflection : public FieldBase
{
public:
    StrattonChuReflection(ISurface& surf, IField& field, const SurfaceRegion& region);

    FieldValue get(const Position& pos) const override;

private:

    VectorComplex subint_E(double x, double y, const Position& r) const;
    VectorComplex subint_B(double x, double y, const Position& r) const;

    ISurface& m_surf;
    IField& m_field;
    SurfaceRegion m_region;
};


#endif // STRATTONCHUFIELD_HPP
