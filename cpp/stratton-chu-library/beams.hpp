#ifndef BEAMS_HPP
#define BEAMS_HPP

#include "types.hpp"

#include <functional>
#include <complex>

constexpr static double c = 29979245800.0;

class ISurface
{
public:
    virtual StaticVector<3> point(const StaticVector<2>& pos) = 0;
    virtual StaticVector<3> dS_over_dxdy(const StaticVector<2>& pos) = 0;
    virtual ~ISurface() = default;
};

struct SurfaceRegion
{
    double x_min = -1.0, x_max = 1.0;
    double y_min = -1.0, y_max = 1.0;
};

struct FieldValue
{
    VectorComplex E;
    VectorComplex B;
};

class IField
{
public:
    virtual FieldValue get(const Position& pos) = 0;
    virtual double k() = 0;
    virtual double omega() = 0;
    virtual double lambda() = 0;

    virtual ~IField() = default;
};

class FieldBase : public IField
{
public:
    FieldBase(double lambda);

    double k() override;
    double omega() override;
    double lambda() override;

protected:
    double m_lambda, m_k, m_omega;
};

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
/*
 * class PlaneSurface(ISurface):
    def __init__(self, r0, alpha, beta):
        """
        Constructor of Plane surface r = r0 + alpha * x + beta * y
        :param r0: any point of surface
        :param alpha: direction vector 1 (3D)
        :param beta: direction vector 2 (3D)
        """
        self.alpha = alpha
        self.beta = beta
        self.dS_over_dxdy_const = np.cross(alpha, beta)
        self.r0 = r0

    # @jit(nopython=True)
    def r(self, xy):
        return self.r0 + self.alpha * xy[0] + self.beta * xy[1]

    def dS_over_dxdy(self, xy):
        return self.dS_over_dxdy_const
        */

class PlaneSurface : public ISurface
{
public:
    PlaneSurface(const StaticVector<3>& r0, const StaticVector<3>& alpha, const StaticVector<3>& beta);
    StaticVector<3> point(const StaticVector<2>& pos) override;
    StaticVector<3> dS_over_dxdy(const StaticVector<2>& pos) override;

private:
    StaticVector<3> r0;
    StaticVector<3> alpha, beta;
    StaticVector<3> dS_over_dxdy_const;
};


#endif // BEAMS_HPP
