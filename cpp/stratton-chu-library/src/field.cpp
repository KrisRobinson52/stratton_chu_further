#include "stratton-chu/field.hpp"

FieldBase::FieldBase(double lambda) :
    m_lambda(lambda)
{
    m_k = 2 * M_PI / lambda;
    m_omega = m_k * c;
}

double FieldBase::k()
{
    return m_k;
}

double FieldBase::omega()
{
    return m_omega;
}

double FieldBase::lambda()
{
    return m_lambda;
}
