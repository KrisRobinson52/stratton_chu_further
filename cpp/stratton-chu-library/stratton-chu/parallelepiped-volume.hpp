#ifndef PARALLELEPIPEDVOLUME_HPP
#define PARALLELEPIPEDVOLUME_HPP

#include "stratton-chu/volume.hpp"

class PrlppdVolume : public IVolume
{
public:
    PrlppdVolume(const Position& r0, const Vector& alpha1, const Vector& alpha2, const Vector& beta);

    Position point(const Vector& pos) const override;
    Vector tau1(const Vector& pos) const override;
    Vector tau2(const Vector& pos) const override;
    Vector tau3(const Vector& pos) const override;

private:
    Position m_r0;
    Vector m_alpha1, m_alpha2, m_beta;
};

#endif // PARALLELEPIPEDVOLUME_HPP
