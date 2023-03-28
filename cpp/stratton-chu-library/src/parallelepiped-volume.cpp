#include "stratton-chu/parallelepiped-volume.hpp"


PrlppdVolume::PrlppdVolume(const Position& r0, const Vector& alpha1, const Vector& alpha2, const Vector& beta) :
    m_r0(r0), m_alpha1(alpha1), m_alpha2(alpha2), m_beta(beta)
{
    
}

Position PrlppdVolume::point(const Vector& pos) const
{
    Position result;
    result = m_r0 + m_alpha1 * pos[0] + m_alpha2 * pos[1] + m_beta * pos[2];
    return result;
}

Vector PrlppdVolume::tau1(const Vector& pos) const
{
    return m_alpha1;
}

Vector PrlppdVolume::tau2(const Vector& pos) const
{
    return m_alpha2;
}

Vector PrlppdVolume::tau3(const Vector& pos) const
{
    return m_beta;
}
