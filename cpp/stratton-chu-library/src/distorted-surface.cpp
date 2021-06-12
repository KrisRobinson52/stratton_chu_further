#include "stratton-chu/distorted-surface.hpp"

SurfaceDistortion::SurfaceDistortion(
        const ISurface& pure_surface,
        const Vector& v,
        const std::vector<DistortionHarmonic>& harmonics) :  // v - unit vector
    m_pure_surface(pure_surface), m_harmonics(harmonics), m_v(v)
{}

Position SurfaceDistortion::point(const Vector2D& pos) const
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

Vector SurfaceDistortion::tau1(const Vector2D& pos) const
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

Vector SurfaceDistortion::tau2(const Vector2D& pos) const
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
