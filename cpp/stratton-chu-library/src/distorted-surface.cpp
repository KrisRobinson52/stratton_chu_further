#include "stratton-chu/distorted-surface.hpp"

#include <boost/math/special_functions/legendre.hpp>

using namespace boost::math;

double legendre(int num, double arg)
{
    if (arg > 1.0)
        return legendre_p(num, 1.0);
    if (arg < -1.0)
        return legendre_p(num, -1.0);
    return legendre_p(num, arg);
}


double legendre_derivative(int num, double arg)
{
    if (arg > 1.0)
        return 0.0;
    if (arg < -1.0)
        return 0.0;
    return legendre_p_prime(num, arg);
}

SurfaceDistortionHarmonic::SurfaceDistortionHarmonic(
        const ISurface& pure_surface,
        const Vector& v,
        const std::vector<DistortionHarmonic>& harmonics) :  // v - unit vector
    m_pure_surface(pure_surface), m_harmonics(harmonics), m_v(v)
{}

Position SurfaceDistortionHarmonic::point(const Vector2D& pos) const
{
    Position result = m_pure_surface.point(pos);
    for (size_t i = 0; i < m_harmonics.size(); i++)
    {
        double shift = m_harmonics[i].ampl * cos(m_harmonics[i].kx * pos[0] + m_harmonics[i].ky * pos[1]);
        Vector delta = m_v * shift;
        result += delta;
    }

    return result;
}

Vector SurfaceDistortionHarmonic::tau1(const Vector2D& pos) const
{
    Vector result = m_pure_surface.tau1(pos);
    for (size_t i = 0; i < m_harmonics.size(); i++)
    {
        double shift = - m_harmonics[i].ampl * m_harmonics[i].kx * sin(m_harmonics[i].kx * pos[0] + m_harmonics[i].ky * pos[1]);
        Vector delta =  m_v * shift;
        result += delta;
    }
    return result;
}

Vector SurfaceDistortionHarmonic::tau2(const Vector2D& pos) const
{
    Vector result = m_pure_surface.tau2(pos);
    for (size_t i = 0; i < m_harmonics.size(); i++)
    {
        double shift = - m_harmonics[i].ampl * m_harmonics[i].ky * sin(m_harmonics[i].kx * pos[0] + m_harmonics[i].ky * pos[1]);
        Vector delta = m_v * shift;
        result += delta;
    }
    return result;
}

SurfaceDistortionLegendre::DistortionPolinom::DistortionPolinom(double ampl, double alpha, int number) :
    ampl(ampl), direction(cos(alpha), sin(alpha)), number(number)
{
}

SurfaceDistortionLegendre::SurfaceDistortionLegendre(
        const ISurface& pure_surface,
        const Vector& v,
        double radius,
        Vector2D center,
        const std::vector<DistortionPolinom>& harmonics) :  // v - unit vector
    m_pure_surface(pure_surface), m_harmonics(harmonics), m_v(v), m_center(center), m_radius(radius)
{
}

Position SurfaceDistortionLegendre::point(const Vector2D& pos) const
{
    Position result = m_pure_surface.point(pos);
    for (size_t i = 0; i < m_harmonics.size(); i++)
    {
        double shift = m_harmonics[i].ampl * legendre(m_harmonics[i].number, m_harmonics[i].direction * (pos - m_center) / m_radius);
        Vector delta = m_v * shift;
        result += delta;
    }

    return result;
}

Vector SurfaceDistortionLegendre::tau1(const Vector2D& pos) const
{
    Vector result = m_pure_surface.tau1(pos);
    for (size_t i = 0; i < m_harmonics.size(); i++)
    {
        double shift = m_harmonics[i].ampl * m_harmonics[i].direction[0] / m_radius * legendre_derivative(m_harmonics[i].number, m_harmonics[i].direction * (pos - m_center) / m_radius);
        Vector delta =  m_v * shift;
        result += delta;
    }
    return result;
}

Vector SurfaceDistortionLegendre::tau2(const Vector2D& pos) const
{
    Vector result = m_pure_surface.tau2(pos);
    for (size_t i = 0; i < m_harmonics.size(); i++)
    {
        double shift = m_harmonics[i].ampl * m_harmonics[i].direction[1] / m_radius * legendre_derivative(m_harmonics[i].number, m_harmonics[i].direction * (pos - m_center) / m_radius);
        Vector delta =  m_v * shift;
        result += delta;
    }
    return result;
}

