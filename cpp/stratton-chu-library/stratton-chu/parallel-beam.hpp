#ifndef PARALLELBEAM_HPP
#define PARALLELBEAM_HPP

#include "stratton-chu/field.hpp"

#include <functional>

class ParallelBeamZ : public FieldBase
{
public:
    using FieldAmpl = std::function<double(double, double)>;

    ParallelBeamZ(double lambda, FieldAmpl Ex, FieldAmpl Ey, double z0);

    FieldValue get(const Position& pos) override;

private:
    FieldAmpl m_Ex, m_Ey;
    double m_z0;
};

#endif // PARALLELBEAM_HPP
