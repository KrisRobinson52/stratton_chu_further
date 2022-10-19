#include "stratton-chu/spec-inverse-fourier.hpp"

#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <array>

// using namespace std::complex_literals;

//template<std::size_t SIZE>

SpecIFT::SpecIFT(std::vector<StrattonChuReflection> &refs, std::vector<double> &freqs, std::vector<std::complex<double>> &amplF) :
    m_refs(refs),
    m_freqs(freqs), m_amplF(amplF)
{
}

FieldValue SpecIFT::get(const Position& pos, const double time) const 
{
    FieldValue result;
    //std::cout << "refs "<<m_refs.size() << std::endl;

    //std::cout << "freqs "<< m_freqs.size() << std::endl;

    //std::cout << "ampls "<< m_amplF.size() << std::endl;
    //std::cout << m_freqs[0]<< std::endl;
    //std::cout << m_amplF[0]<< std::endl;

    for (size_t i=0; i<m_refs.size(); i++)
    {   
        //m_refs[i].get(pos);
        //std::cout << i << std::endl;
        result.E=result.E+m_refs[i].get(pos).E*m_amplF[i]*exp(std::complex<double>(0.0, 1.0) * m_freqs[i] * time);
    }
    //std::cout << pos[0] << " " << pos[1]<< " "<< pos[2]<< std::endl;
    return result;
}