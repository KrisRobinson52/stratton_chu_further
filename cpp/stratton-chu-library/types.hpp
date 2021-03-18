#ifndef TYPES_HPP
#define TYPES_HPP

#include "geometry.hpp"

#include <complex>

using Complex = std::complex<double>;
using Position = StaticVector<3>;
using Vector = StaticVector<3>;
using Vector2D = StaticVector<2>;
using VectorComplex = StaticVector<3, Complex>;

#endif // TYPES_HPP
