#include "stratton-chu/field.hpp"

FieldBase::FieldBase(double lambda) :
    m_lambda(lambda)
{
    m_k = 2 * M_PI / lambda;
    m_omega = m_k * c;
}

double FieldBase::k() const
{
    return m_k;
}

double FieldBase::omega() const
{
    return m_omega;
}

double FieldBase::lambda() const
{
    return m_lambda;
}

SummaryField::SummaryField(const IField& field1, const IField& field2) :
    FieldBase(field1.lambda()), m_field1(field1), m_field2(field2)
{
    if (field1.lambda() != field2.lambda())
        throw std::invalid_argument("Cannot sum fields with different wave length");
}

FieldValue SummaryField::get(const Position& pos) const
{
    return m_field1.get(pos) + m_field2.get(pos);
}

SummaryField operator+(const IField& left, const IField& right)
{
    return SummaryField(left, right);
}
