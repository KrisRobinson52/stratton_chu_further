#ifndef PARALLELBEAM_HPP
#define PARALLELBEAM_HPP

#include "stratton-chu/field.hpp"

#include <functional>

using FieldAmpl = std::function<Complex(double, double)>;

class ParallelBeamZ : public FieldBase  // Тут нет расплывания получается в фазовом множителе...
{
public:
    ParallelBeamZ(double lambda, FieldAmpl Ex, FieldAmpl Ey, double z0);

    FieldValue get(const Position& pos) const override;

private:
    FieldAmpl m_Ex, m_Ey;
    double m_z0;
};

class ParallelBeamAlpha : public FieldBase
{
public:
    /* alpha1 & alpha2 define beam plane
     * [alpha1 x alpha2] defines beam direction
     *
     * E1 - complex amplitude of the field polarized along alpha1
     * E2 - complex amplitude of the field polarized along alpha2
     * E1 and E2 expected to be functions defined in the beam plane
     */
    ParallelBeamAlpha(double lambda, const Position& r0, const Vector& alpha1, const Vector& alpha2,
                      const FieldAmpl& E1, const FieldAmpl& E2);

    FieldValue get(const Position& pos) const override;

private:
    Position m_r0;
    Vector m_alpha1, m_alpha2;
    Vector m_alpha;
    FieldAmpl m_E1, m_E2;
};

#endif // PARALLELBEAM_HPP
