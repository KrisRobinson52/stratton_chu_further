#ifndef SPECINVERSEFOURIER_HPP
#define SPECINVERSEFOURIER_HPP

#include "stratton-chu/types.hpp"
#include "stratton-chu/field.hpp"
#include <vector>
#include <complex>
#include "stratton-chu/stratton-chu-field.hpp"

//template<std::size_t SIZE>

class SpecIFT
{
public:
    SpecIFT(std::vector<StrattonChuReflection> &refs, std::vector<double> &freqs, std::vector<std::complex<double>> &amplF);

    FieldValue get(const Position& pos, const double time) const;

private:
    std::vector<StrattonChuReflection> &m_refs; 
    std::vector<double> &m_freqs;
    std::vector<std::complex<double>> &m_amplF;
};


#endif //SPECINVERSEFOURIER_HPP
