#ifndef FIELDINREGION_HPP
#define FIELDINREGION_HPP

#include "stratton-chu/types.hpp"
#include "stratton-chu/field.hpp"
#include "stratton-chu/surface.hpp"
#include <iostream>
#include "stratton-chu/stratton-chu-field.hpp"
#include "stratton-chu/spec-inverse-fourier.hpp"

class FieldInRegion
{
public:
    FieldInRegion(ISurface& surf, SurfaceRegion& region, std::vector<StrattonChuReflection> &fields, int n_harmonics, int n_points);

    //std::vector<std::vector<FieldValue>> InRegion(const ISurface& surface, const SurfaceRegion& region, int n_points);

    FieldValue get(int ifreq, int irow, int icolumn) const;
    std::vector<std::vector<FieldValue>> getIFT(std::vector<double> &freqs, std::vector<std::complex<double>> &amplF, const double time);

private:

    std::vector<std::vector<std::vector<FieldValue>>> &m_FieldCatalogue;
    ISurface& m_surf;
    std::vector<StrattonChuReflection> &m_fields; 
    SurfaceRegion m_region;
    int mn_harm;
    int mn_points;
};


#endif // FIELDINREGION_HPP