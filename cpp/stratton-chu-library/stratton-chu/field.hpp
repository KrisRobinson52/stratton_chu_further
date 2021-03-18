#ifndef FIELD_HPP
#define FIELD_HPP

#include "stratton-chu/types.hpp"

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

#endif // FIELD_HPP
