#ifndef STRATTONCHUFIELD_HPP
#define STRATTONCHUFIELD_HPP

#include "stratton-chu/types.hpp"
#include "stratton-chu/field.hpp"
#include "stratton-chu/surface.hpp"
#include "stratton-chu/volume.hpp"
#include <iostream>

class StrattonChuReflection : public FieldBase
{
public:
    StrattonChuReflection(ISurface& surf, IField& field, const SurfaceRegion& region);
    // StrattonChuReflection(ISurface& surf, auto& fields, const SurfaceRegion& region);
    // StrattonChuReflection(IVolume& vol, auto& fields, const VolumeRegion& region3d);

    //std::vector<std::vector<FieldValue>> InRegion(const ISurface& surface, const SurfaceRegion& region, int n_points);

    FieldValue get(const Position& pos) const override;

private:

    VectorComplex subint_E(double x, double y, const Position& r) const;
    VectorComplex subint_B(double x, double y, const Position& r) const;

    ISurface& m_surf;
    IField& m_field;
    SurfaceRegion m_region;

    // IVolume& m_vol;
    // VolumeRegion& m_region3d;

};

//для тестов
class StrattonChuReflectionLowPrec : public FieldBase
{
public:
    StrattonChuReflectionLowPrec(ISurface& surf, IField& field, const SurfaceRegion& region);
    // StrattonChuReflection(ISurface& surf, auto& fields, const SurfaceRegion& region);
    // StrattonChuReflection(IVolume& vol, auto& fields, const VolumeRegion& region3d);

    //std::vector<std::vector<FieldValue>> InRegion(const ISurface& surface, const SurfaceRegion& region, int n_points);

    FieldValue get(const Position& pos) const override;

private:

    VectorComplex subint_E(double x, double y, const Position& r) const;
    VectorComplex subint_EE(double x, double y, const Position& r, int& checker) const;
    VectorComplex subint_B(double x, double y, const Position& r) const;
    VectorComplex subint_BB(double x, double y, const Position& r, int& checker) const;

    ISurface& m_surf;
    IField& m_field;
    SurfaceRegion m_region;

    // IVolume& m_vol;
    // VolumeRegion& m_region3d;

};


#endif // STRATTONCHUFIELD_HPP
