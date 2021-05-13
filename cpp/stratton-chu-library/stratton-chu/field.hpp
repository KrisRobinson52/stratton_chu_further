#ifndef FIELD_HPP
#define FIELD_HPP

#include "stratton-chu/types.hpp"

constexpr double c = 29979245800.0;

struct FieldValue
{
    VectorComplex E;
    VectorComplex B;

    FieldValue operator+(const FieldValue& right)
    {
        FieldValue result(*this);
        result.E += right.E;
        result.B += right.B;
        return result;
    }
};

class IField
{
public:
    virtual FieldValue get(const Position& pos) const = 0;
    virtual double k() const = 0;
    virtual double omega() const = 0;
    virtual double lambda() const = 0;

    virtual ~IField() = default;
};

class FieldBase : public IField
{
public:
    FieldBase(double lambda);

    double k() const override;
    double omega() const override;
    double lambda() const override;

protected:
    double m_lambda, m_k, m_omega;
};


class SummaryField : public FieldBase
{
public:
    SummaryField(const IField& field1, const IField& field2);
    FieldValue get(const Position& pos) const override;

private:
    const IField& m_field1;
    const IField& m_field2;
};

SummaryField operator+(const IField& left, const IField& right);

#endif // FIELD_HPP
