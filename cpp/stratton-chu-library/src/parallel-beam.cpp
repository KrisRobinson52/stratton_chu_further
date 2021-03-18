#include "stratton-chu/parallel-beam.hpp"

using namespace std::complex_literals;

ParallelBeamZ::ParallelBeamZ(double lambda, FieldAmpl Ex, FieldAmpl Ey, double z0) :
    FieldBase(lambda), m_Ex(Ex), m_Ey(Ey), m_z0(z0)
{
}

FieldValue ParallelBeamZ::get(const Position& pos)
{
    double x = pos[0], y = pos[1], z = pos[2];
    std::complex phase_mul = std::exp( -1i * m_k *(z - m_z0));
    FieldValue result;

    std::complex Ex = m_Ex(x, y) * phase_mul;
    std::complex Ey = m_Ey(x, y) * phase_mul;

    result.E[0] = Ex; result.E[0] = Ey; result.E[0] = 0;
    result.B[0] = -Ey; result.B[0] = Ex; result.E[0] = 0;

    return result;
}
