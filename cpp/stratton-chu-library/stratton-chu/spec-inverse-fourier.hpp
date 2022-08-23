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
    //SpecIFT(std::vector<StrattonChuReflection> &refs, std::array <double, SIZE> &freqs, std::array <std::complex<double>, SIZE> &amplF);

    FieldValue get(const Position& pos, const double time) const;

private:
    std::vector<StrattonChuReflection> &m_refs; 
    std::vector<double> &m_freqs;
    //std::array <double, harmonics_number> &m_freqs;
    std::vector<std::complex<double>> &m_amplF;
    //std::array <std::complex<double>, harmonics_number> &m_amplF;
};


#endif //SPECINVERSEFOURIER_HPP
