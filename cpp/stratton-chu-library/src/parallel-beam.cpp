#include "stratton-chu/parallel-beam.hpp"

#include <cmath>

// using namespace std::complex_literals;

ParallelBeamZ::ParallelBeamZ(double lambda, FieldAmpl Ex, FieldAmpl Ey, double z0) :
    FieldBase(lambda), m_Ex(Ex), m_Ey(Ey), m_z0(z0)
{
}

FieldValue ParallelBeamZ::get(const Position& pos)
{
    double x = pos[0], y = pos[1], z = pos[2];
    Complex phase_mul = std::exp( - Complex(0.0, 1.0) * m_k *(z - m_z0));

    FieldValue result;

    std::complex Ex = m_Ex(x, y) * phase_mul;
    std::complex Ey = m_Ey(x, y) * phase_mul;

    result.E[0] = Ex; result.E[1] = Ey; result.E[2] = 0;
    result.B[0] = Ey; result.B[1] = - Ex; result.B[2] = 0;

    return result;
}

ParallelBeamAlpha::ParallelBeamAlpha(double lambda, const Position& r0, const Vector& alpha1, const Vector& alpha2,
                                     const FieldAmpl& E1, const FieldAmpl& E2) :
    FieldBase(lambda), m_r0(r0),
    m_alpha1(alpha1 / alpha1.norm()), m_alpha2(alpha2 / alpha2.norm()),
    m_alpha(m_alpha1 % m_alpha2), m_E1(E1), m_E2(E2)
{
}

FieldValue ParallelBeamAlpha::get(const Position& pos)
{
    Vector delta_r = (pos - m_r0);
    double x1 = delta_r * m_alpha1,
           x2 = delta_r * m_alpha2,
           l = delta_r * m_alpha ;
    Complex phase_mul = std::exp( Complex(0.0, 1.0) * m_k * l );

    FieldValue result;
    Complex E1 = m_E1(x1, x2) * phase_mul;
    Complex E2 = m_E2(x1, x2) * phase_mul;

    result.E = E1 * m_alpha1.cast<Complex>() + E2 * m_alpha2.cast<Complex>();
    result.B = E1 * m_alpha2.cast<Complex>() - E2 * m_alpha1.cast<Complex>();

    return result;
}
