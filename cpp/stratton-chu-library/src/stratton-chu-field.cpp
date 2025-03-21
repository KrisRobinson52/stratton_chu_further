#include "stratton-chu/stratton-chu-field.hpp"

#include "stratton-chu/utils.hpp"
#include <iostream>


#include <functional>
#include <unordered_map>
#include "constants.hpp"
#include "cache.hpp"
#include "static_local_tracker.hpp"

#define GET_VARIABLE_NAME(var) (#var)


using namespace std::complex_literals;
using namespace Caching;


StrattonChuReflection::StrattonChuReflection(ISurface& surf, IField& field, const SurfaceRegion& region) :
    FieldBase(field.lambda()),
    m_surf(surf), m_field(field), m_region(region)
{
}

// StrattonChuReflection::StrattonChuReflection(ISurface& surf, auto& fields, const SurfaceRegion& region) :
//     FieldBase(fields.getLambda()),
//     m_surf(surf), m_region(region)
// {
// }

// StrattonChuReflection::StrattonChuReflection(IVolume& vol, auto& fields, const VolumeRegion& region3d) :
//     FieldBase(fields.getLambda()),
//     m_vol(vol), m_region3d(region3d)
// {
// }

FieldValue StrattonChuReflection::get(const Position& pos) const
{
    //static std::unordered_map<Key, FieldValue> cache{(focal_points*focal_points+focal_points_transversal*(focal_points_transversal-1))*harmonics_number};
    
    // Dependances deps{pos[0], pos[1], pos[2], this->lambda()};

    FieldValue result;
    // static ConcurrentCache<decltype(deps), FieldValue> conc_cache{};
    // STATIC_LOCAL_TRACK(conc_cache);
    // auto stored = conc_cache.load(deps);
/*
    result.E = integrate_trapezoid<VectorComplex>(
        [this, &pos] (double x, double y) -> VectorComplex { return subint_E(x, y, pos); },
        m_region,
        200, 200
    );
*/
    // std::cout << "in" << std::endl;

    // if (stored.has_value())
    // {        
    //     //[[maybe_unused]] static const auto flag = [&]{conc_cache.set_stores_availability(false); return true;}();
    //     result = stored.value();
    // }
    // else
    // {

        result.E = integrate_cubature(
            [this, &pos] (double x, double y) -> VectorComplex { return subint_E(x, y, pos); },
            m_region, 1e-3, 1e-3);
            // m_region, 1e-1, 1e-1);
        result.E *= 1 / (4 * M_PI);


    // std::cout << "out" << std::endl;
    //     conc_cache.store(deps, result);
    // }

        result.B = integrate_cubature(
            [this, &pos] (double x, double y) -> VectorComplex { return subint_B(x, y, pos); },
            m_region, 1e-3, 1e-3);
            // m_region, 1e-1, 1e-1);
        result.B *= 1 / (4 * M_PI);

    // @TODO Calculate B
    //cache[key] = result;
    return result;
}

VectorComplex StrattonChuReflection::subint_E(double x, double y, const Position& r) const
{
    Vector2D xy = {x, y};
    Vector r0 = m_surf.point(xy);
    VectorComplex N = m_surf.dS_over_dxdy(xy).cast<Complex>();

    // std::cout << "zdes1" << std::endl;

    FieldValue field = m_field.get(r0);

    // std::cout << "zdes12" << std::endl;

    GreenFunc green(r, r0, m_k);

    //std::cout << "zdes2" << std::endl;

    VectorComplex first_term = (N % field.B);
    first_term *= 2.0 * Complex(0.0, 1.0) * m_k * green.value();

    // std::cout << "zdes3" << std::endl;

    VectorComplex second_term = green.gradient();
    second_term *= (2.0 * (N.operator*(field.E) ));

    // std::cout << "4 " << std::endl;
    // std::cout << second_term[0] << std::endl;

    return first_term + second_term;
}


VectorComplex StrattonChuReflection::subint_B(double x, double y, const Position& r) const
{
    Vector2D xy = {x, y};
    Vector r0 = m_surf.point(xy);
    VectorComplex N = m_surf.dS_over_dxdy(xy).cast<Complex>();

    FieldValue field = m_field.get(r0);
    GreenFunc green(r, r0, m_k);

    VectorComplex mul1 = (N % field.B);
    VectorComplex mul2 = green.gradient();
    VectorComplex mul = (mul1 % mul2);
    mul *= 2.0;

    return mul;
    
    // return VectorComplex();
}



StrattonChuReflectionLowPrec::StrattonChuReflectionLowPrec(ISurface& surf, IField& field, const SurfaceRegion& region) :
    FieldBase(field.lambda()),
    m_surf(surf), m_field(field), m_region(region)
{
}

FieldValue StrattonChuReflectionLowPrec::get(const Position& pos) const
{
    //static std::unordered_map<Key, FieldValue> cache{(focal_points*focal_points+focal_points_transversal*(focal_points_transversal-1))*harmonics_number};
    
    // Dependances deps{pos[0], pos[1], pos[2], this->lambda()};

    FieldValue result;
    // static ConcurrentCache<decltype(deps), FieldValue> conc_cache{};
    // STATIC_LOCAL_TRACK(conc_cache);
    // auto stored = conc_cache.load(deps);
/*
    result.E = integrate_trapezoid<VectorComplex>(
        [this, &pos] (double x, double y) -> VectorComplex { return subint_E(x, y, pos); },
        m_region,
        200, 200
    );
*/
    std::cout << "in" << std::endl;

    // if (stored.has_value())
    // {        
    //     //[[maybe_unused]] static const auto flag = [&]{conc_cache.set_stores_availability(false); return true;}();
    //     result = stored.value();
    // }
    // else
    // {

    int checkerE=0, checkerB=0;

        // result.E = integrate_cubature(
        //     [this, &pos, &checkerE] (double x, double y) -> VectorComplex { return subint_EE(x, y, pos, checkerE); },
        //     m_region, 1e-3, 1e-3);
        //     // m_region, 1e-1, 1e-1);
        // result.E *= 1 / (4 * M_PI);

        result.E = integrate_cubature_low_prec(
            [this, &pos] (double x, double y) -> VectorComplex { /*std::cout << "in E" << std::endl;*/ return subint_E(x, y, pos); },
            // m_region, 1e-3, 1e-3);
            // m_region, 1e-2, 1e-2);
            m_region, 1e-1, 1e-1);
        result.E *= 1 / (4 * M_PI);

// 
    // std::cout << "outE" << std::endl;
    //     conc_cache.store(deps, result);
    // }

        // result.B = integrate_cubature_low_prec(
        //     [this, &pos] (double x, double y) -> VectorComplex { return subint_B(x, y, pos); },
        //     // m_region, 1e-3, 1e-3);
        //     m_region, 1e-1, 1e-1);
        // result.B *= 1 / (4 * M_PI);

        // result.B = integrate_cubature_low_prec(
        //     [this, &pos, &checkerB] (double x, double y) -> VectorComplex { return subint_BB(x, y, pos, checkerB); },
        //     m_region, 1e-3, 1e-3);
        //     // m_region, 1e-1, 1e-1);
        // result.B *= 1 / (4 * M_PI);

    // @TODO Calculate B
    //cache[key] = result;

    // std::cout << "strat points E: " << checkerE  << ", B: " << checkerB << std::endl;
    return result;
}

VectorComplex StrattonChuReflectionLowPrec::subint_E(double x, double y, const Position& r) const
{
    // std::cout << "in in " << std::endl;
    Vector2D xy = {x, y};
    Vector r0 = m_surf.point(xy);
    VectorComplex N = m_surf.dS_over_dxdy(xy).cast<Complex>();

    // std::cout << "zdes1" << std::endl;

    FieldValue field = m_field.get(r0);


    // std::cout << "zdes12" << std::endl;

    GreenFunc green(r, r0, m_k);

    //std::cout << "zdes2" << std::endl;

    VectorComplex first_term = (N % field.B);
    first_term *= 2.0 * Complex(0.0, 1.0) * m_k * green.value();

    // std::cout << "zdes3" << std::endl;

    VectorComplex second_term = green.gradient();
    second_term *= (2.0 * (N.operator*(field.E) ));

    // std::cout << "4 " << std::endl;
    // std::cout << second_term[0] << std::endl;

    return first_term + second_term;
}

VectorComplex StrattonChuReflectionLowPrec::subint_EE(double x, double y, const Position& r, int& checker) const
{
    Vector2D xy = {x, y};
    Vector r0 = m_surf.point(xy);
    VectorComplex N = m_surf.dS_over_dxdy(xy).cast<Complex>();

    // std::cout << "zdes1" << std::endl;

    FieldValue field = m_field.get(r0);

    // std::cout << "zdes12" << std::endl;

    GreenFunc green(r, r0, m_k);

    //std::cout << "zdes2" << std::endl;

    VectorComplex first_term = (N % field.B);
    first_term *= 2.0 * Complex(0.0, 1.0) * m_k * green.value();

    // std::cout << "zdes3" << std::endl;

    VectorComplex second_term = green.gradient();
    second_term *= (2.0 * (N.operator*(field.E) ));

    // std::cout << "4 " << std::endl;
    // std::cout << second_term[0] << std::endl;
    checker++;

    return first_term + second_term;
}

VectorComplex StrattonChuReflectionLowPrec::subint_B(double x, double y, const Position& r) const
{
    Vector2D xy = {x, y};
    Vector r0 = m_surf.point(xy);
    VectorComplex N = m_surf.dS_over_dxdy(xy).cast<Complex>();

    FieldValue field = m_field.get(r0);

    GreenFunc green(r, r0, m_k);

    VectorComplex mul1 = (N % field.B);
    VectorComplex mul2 = green.gradient();
    VectorComplex mul = (mul1 % mul2);
    mul *= 2.0;

    return mul;
    
    // return VectorComplex();
}

VectorComplex StrattonChuReflectionLowPrec::subint_BB(double x, double y, const Position& r, int& checker) const
{
    Vector2D xy = {x, y};
    Vector r0 = m_surf.point(xy);
    VectorComplex N = m_surf.dS_over_dxdy(xy).cast<Complex>();

    FieldValue field = m_field.get(r0);

    GreenFunc green(r, r0, m_k);

    VectorComplex mul1 = (N % field.B);
    VectorComplex mul2 = green.gradient();
    VectorComplex mul = (mul1 % mul2);
    mul *= 2.0;

    checker++;

    return mul;
    
    // return VectorComplex();
}
