#ifndef FIELD_HPP
#define FIELD_HPP

#include "stratton-chu/types.hpp"
#include <bit>
#include <array>
#include <span>

#include "stdafx.h"
#include "interpolation.h"
#include "stratton-chu/parabolic-surface.hpp"
#include "stratton-chu/elliptic-surface.hpp"
#include "threading.hpp"
// #include "constants.hpp"

constexpr double c = 29979245800.0;

struct FieldValue
{
    VectorComplex E;
    VectorComplex B;

    auto serialize () const {
        std::array<std::byte, sizeof(FieldValue)> bytes;
        std::memcpy(bytes.data(),this,sizeof(FieldValue));
        return bytes;
    }

    static auto deserialize (std::span<std::byte> bytes) {
        FieldValue fieldv;
        std::memcpy(&fieldv,bytes.data(),sizeof(FieldValue));
        return fieldv;
    }


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

template <const size_t n_points>
class FieldInterpolated : public FieldBase
{
public:
    FieldInterpolated(double lambda, std::array<double, n_points> &xs, std::array<double, n_points> &ys, std::array<double, n_points> &zs, double PointsDist, std::array<FieldValue, n_points> &FieldInRegion) :
        FieldBase(lambda), m_npoints(n_points), m_xs(xs),  m_ys(ys),  m_zs(zs), m_FieldInRegion(FieldInRegion), m_PointsDist(PointsDist)
    {
        // alglib::rbfcreate(3,3, m_EFieldReal);
        // alglib::rbfcreate(3,3, m_EFieldImag);
        // alglib::rbfcreate(3,3, m_BFieldReal);
        // alglib::rbfcreate(3,3, m_BFieldImag);

        alglib::rbfcreate(3,3, m_EFieldMod);
        alglib::rbfcreate(3,3, m_EFieldPhase);
        alglib::rbfcreate(3,3, m_BFieldMod);
        alglib::rbfcreate(3,3, m_BFieldPhase);

        // m_Points_and_EReals.setlength(m_npoints, 6);
        // m_Points_and_EImags.setlength(m_npoints, 6);
        // m_Points_and_BReals.setlength(m_npoints, 6);
        // m_Points_and_BImags.setlength(m_npoints, 6);

        m_Points_and_EMods.setlength(m_npoints, 6);
        m_Points_and_EPhases.setlength(m_npoints, 6);
        m_Points_and_BMods.setlength(m_npoints, 6);
        m_Points_and_BPhases.setlength(m_npoints, 6);
       
        // for (int i=0; i<m_npoints; i++) {
        //     m_Points_and_EReals[i][0]=m_xs[i];
        //     m_Points_and_EReals[i][1]=m_ys[i];
        //     m_Points_and_EReals[i][2]=m_zs[i];
        //     m_Points_and_EReals[i][3]=vec_real(m_FieldInRegion[i].E)[0];
        //     m_Points_and_EReals[i][4]=vec_real(m_FieldInRegion[i].E)[1];
        //     m_Points_and_EReals[i][5]=vec_real(m_FieldInRegion[i].E)[2];

        //     m_Points_and_EImags[i][0]=m_xs[i];
        //     m_Points_and_EImags[i][1]=m_ys[i];
        //     m_Points_and_EImags[i][2]=m_zs[i];
        //     m_Points_and_EImags[i][3]=vec_imag(m_FieldInRegion[i].E)[0];
        //     m_Points_and_EImags[i][4]=vec_imag(m_FieldInRegion[i].E)[1];
        //     m_Points_and_EImags[i][5]=vec_imag(m_FieldInRegion[i].E)[2];

        //     m_Points_and_BReals[i][0]=m_xs[i];
        //     m_Points_and_BReals[i][1]=m_ys[i];
        //     m_Points_and_BReals[i][2]=m_zs[i];
        //     m_Points_and_BReals[i][3]=vec_real(m_FieldInRegion[i].B)[0];
        //     m_Points_and_BReals[i][4]=vec_real(m_FieldInRegion[i].B)[1];
        //     m_Points_and_BReals[i][5]=vec_real(m_FieldInRegion[i].B)[2];

        //     m_Points_and_BImags[i][0]=m_xs[i];
        //     m_Points_and_BImags[i][1]=m_ys[i];
        //     m_Points_and_BImags[i][2]=m_zs[i];
        //     m_Points_and_BImags[i][3]=vec_imag(m_FieldInRegion[i].B)[0];
        //     m_Points_and_BImags[i][4]=vec_imag(m_FieldInRegion[i].B)[1];
        //     m_Points_and_BImags[i][5]=vec_imag(m_FieldInRegion[i].B)[2];
        // }

        for (int i=0; i<m_npoints; i++) {
            m_Points_and_EMods[i][0]=m_xs[i];
            m_Points_and_EMods[i][1]=m_ys[i];
            m_Points_and_EMods[i][2]=m_zs[i];
            m_Points_and_EMods[i][3]=vec_modulus(m_FieldInRegion[i].E)[0];
            m_Points_and_EMods[i][4]=vec_modulus(m_FieldInRegion[i].E)[1];
            m_Points_and_EMods[i][5]=vec_modulus(m_FieldInRegion[i].E)[2];

            m_Points_and_EPhases[i][0]=m_xs[i];
            m_Points_and_EPhases[i][1]=m_ys[i];
            m_Points_and_EPhases[i][2]=m_zs[i];
            m_Points_and_EPhases[i][3]=vec_phases(m_FieldInRegion[i].E)[0];
            m_Points_and_EPhases[i][4]=vec_phases(m_FieldInRegion[i].E)[1];
            m_Points_and_EPhases[i][5]=vec_phases(m_FieldInRegion[i].E)[2];

            m_Points_and_BMods[i][0]=m_xs[i];
            m_Points_and_BMods[i][1]=m_ys[i];
            m_Points_and_BMods[i][2]=m_zs[i];
            m_Points_and_BMods[i][3]=vec_modulus(m_FieldInRegion[i].B)[0];
            m_Points_and_BMods[i][4]=vec_modulus(m_FieldInRegion[i].B)[1];
            m_Points_and_BMods[i][5]=vec_modulus(m_FieldInRegion[i].B)[2];

            m_Points_and_BPhases[i][0]=m_xs[i];
            m_Points_and_BPhases[i][1]=m_ys[i];
            m_Points_and_BPhases[i][2]=m_zs[i];
            m_Points_and_BPhases[i][3]=vec_phases(m_FieldInRegion[i].B)[0];
            m_Points_and_BPhases[i][4]=vec_phases(m_FieldInRegion[i].B)[1];
            m_Points_and_BPhases[i][5]=vec_phases(m_FieldInRegion[i].B)[2];
        }
        // alglib::rbfsetpoints(m_EFieldReal, m_Points_and_EReals);
        // alglib::rbfsetpoints(m_EFieldImag, m_Points_and_EImags);
        // alglib::rbfsetpoints(m_BFieldReal, m_Points_and_BReals);
        // alglib::rbfsetpoints(m_BFieldImag, m_Points_and_BImags);

        // alglib::rbfsetalgohierarchical(m_EFieldReal, m_PointsDist/4, 2, 0.0);
        // // alglib::rbfsetalgobiharmonic(m_EFieldReal, 1.0E-3);
        // // alglib::rbfsetalgomultiquadricauto(m_EFieldReal, 1.0E-3);
        // alglib::rbfbuildmodel(m_EFieldReal, m_EFieldRealRep);

        // alglib::rbfsetalgohierarchical(m_EFieldImag, m_PointsDist/4, 2, 0.0);
        // // alglib::rbfsetalgobiharmonic(m_EFieldImag, 1.0E-3);
        // // alglib::rbfsetalgomultiquadricauto(m_EFieldImag, 1.0E-3);
        // alglib::rbfbuildmodel(m_EFieldImag, m_EFieldImagRep);

        // alglib::rbfsetalgohierarchical(m_BFieldReal, m_PointsDist/4, 2, 0.0);
        // // alglib::rbfsetalgobiharmonic(m_BFieldReal, 1.0E-3);
        // // alglib::rbfsetalgomultiquadricauto(m_BFieldReal, 1.0E-3);
        // alglib::rbfbuildmodel(m_BFieldReal, m_BFieldRealRep);

        // alglib::rbfsetalgohierarchical(m_BFieldImag, m_PointsDist/4, 2, 0.0);
        // // alglib::rbfsetalgobiharmonic(m_BFieldImag, 1.0E-3);
        // // alglib::rbfsetalgomultiquadricauto(m_BFieldImag, 1.0E-3);
        // alglib::rbfbuildmodel(m_BFieldImag, m_BFieldImagRep);
        

        // std::cout << "EReal check (should be 1) ---> " << int(m_EFieldRealRep.terminationtype) << std::endl;
        // std::cout << "EImag check (should be 1) ---> " << int(m_EFieldImagRep.terminationtype) << std::endl;
        // std::cout << "BReal check (should be 1) ---> " << int(m_BFieldRealRep.terminationtype) << std::endl;
        // std::cout << "BImag check (should be 1) ---> " << int(m_BFieldImagRep.terminationtype) << std::endl;


        alglib::rbfsetpoints(m_EFieldMod, m_Points_and_EMods);
        alglib::rbfsetpoints(m_EFieldPhase, m_Points_and_EPhases);
        alglib::rbfsetpoints(m_BFieldMod, m_Points_and_BMods);
        alglib::rbfsetpoints(m_BFieldPhase, m_Points_and_BPhases);

        // struct Task{
        //     alglib::rbfmodel &m_EFieldMod, &m_EFieldPhase, &m_BFieldMod, &m_BFieldPhase;
        //     alglib::rbfreport &m_EFieldModRep, &m_EFieldPhaseRep, &m_BFieldModRep, &m_BFieldPhaseRep;
        //     double m_PointsDist;
        //     size_t j;
        //     void operator()(){
        //         if (j==0) {
        //             // alglib::rbfsetalgohierarchical(m_EFieldMod, m_PointsDist, 7, 0.0);
        //             alglib::rbfsetalgohierarchical(m_EFieldMod, m_PointsDist, 10, 0.0);
        //             // alglib::rbfsetalgobiharmonic(m_EFieldMod, 0.0);
        //             // alglib::rbfsetalgomultiquadricauto(m_EFieldMod, 0.0);
        //             // alglib::rbfsetalgothinplatespline(m_EFieldMod, 0.0);
        //             alglib::rbfbuildmodel(m_EFieldMod, m_EFieldModRep);
        //             std::cout << "EMod check (should be 1) ---> " << int(m_EFieldModRep.terminationtype) << std::endl;
        //         };

        //         if (j==1) {
        //             // alglib::rbfsetalgohierarchical(m_EFieldPhase, m_PointsDist, 7, 0.0);
        //             alglib::rbfsetalgohierarchical(m_EFieldPhase, m_PointsDist, 10, 0.0);
        //             // alglib::rbfsetalgobiharmonic(m_EFieldPhase, 0.0);
        //             // alglib::rbfsetalgomultiquadricauto(m_EFieldPhase,0.0);
        //             // alglib::rbfsetalgothinplatespline(m_EFieldPhase,0.0);
        //             alglib::rbfbuildmodel(m_EFieldPhase, m_EFieldPhaseRep);
        //             std::cout << "EPhase check (should be 1) ---> " << int(m_EFieldPhaseRep.terminationtype) << std::endl;
        //         };

        //         if (j==2) {
        //             // alglib::rbfsetalgohierarchical(m_BFieldMod, m_PointsDist, 7, 0.0);
        //             alglib::rbfsetalgohierarchical(m_BFieldMod, m_PointsDist, 10, 0.0);
        //             // alglib::rbfsetalgobiharmonic(m_BFieldMod, 0.0);
        //             // alglib::rbfsetalgomultiquadricauto(m_BFieldMod, 0.0);
        //             // alglib::rbfsetalgothinplatespline(m_BFieldMod, 0.0);
        //             alglib::rbfbuildmodel(m_BFieldMod, m_BFieldModRep);
        //             std::cout << "BMod check (should be 1) ---> " << int(m_BFieldModRep.terminationtype) << std::endl;
        //         };

        //         if (j==3) {
        //             // alglib::rbfsetalgohierarchical(m_BFieldPhase, m_PointsDist, 7, 0.0);
        //             alglib::rbfsetalgohierarchical(m_BFieldPhase, m_PointsDist, 10, 0.0);
        //             // alglib::rbfsetalgobiharmonic(m_BFieldPhase, 0.0);
        //             // alglib::rbfsetalgomultiquadricauto(m_BFieldPhase, 0.0);
        //             // alglib::rbfsetalgothinplatespline(m_BFieldPhase, 0.0);
        //             alglib::rbfbuildmodel(m_BFieldPhase, m_BFieldPhaseRep);
        //             std::cout << "BPhase check (should be 1) ---> " << int(m_BFieldPhaseRep.terminationtype) << std::endl;
        //         };
                
        //     }
        //     operator bool()const{
        //         return true;
        //     }
        // };
        // using Threads = Threading::ThreadPool<ConcurrentQueue<Task>,14>;
        // typename Threads::TaskQueue tasks;

        // for (int j = 0; j < 4; j++)
        // {
        //     tasks.push(Task{m_EFieldMod, m_EFieldPhase, m_BFieldMod, m_BFieldPhase, m_EFieldModRep, m_EFieldPhaseRep, m_BFieldModRep, m_BFieldPhaseRep, m_PointsDist, j});
        // }
        // Threads threads{tasks};
        // threads.run();

        // alglib::rbfsetalgohierarchical(m_EFieldMod, m_PointsDist, 7, 0.0);
        alglib::rbfsetalgohierarchical(m_EFieldMod, m_PointsDist, 5, 0.0);
        // alglib::rbfsetalgobiharmonic(m_EFieldMod, 0.0);
        // alglib::rbfsetalgomultiquadricauto(m_EFieldMod, 0.0);
        // alglib::rbfsetalgothinplatespline(m_EFieldMod, 0.0);
        alglib::rbfbuildmodel(m_EFieldMod, m_EFieldModRep);

        // alglib::rbfsetalgohierarchical(m_EFieldPhase, m_PointsDist, 7, 0.0);
        alglib::rbfsetalgohierarchical(m_EFieldPhase, m_PointsDist/4, 3, 0.0);
        // alglib::rbfsetalgobiharmonic(m_EFieldPhase, 0.0);
        // alglib::rbfsetalgomultiquadricauto(m_EFieldPhase,0.0);
        // alglib::rbfsetalgothinplatespline(m_EFieldPhase,0.0);
        alglib::rbfbuildmodel(m_EFieldPhase, m_EFieldPhaseRep);

        // alglib::rbfsetalgohierarchical(m_BFieldMod, m_PointsDist, 7, 0.0);
        alglib::rbfsetalgohierarchical(m_BFieldMod, m_PointsDist, 5, 0.0);
        // alglib::rbfsetalgobiharmonic(m_BFieldMod, 0.0);
        // alglib::rbfsetalgomultiquadricauto(m_BFieldMod, 0.0);
        // alglib::rbfsetalgothinplatespline(m_BFieldMod, 0.0);
        alglib::rbfbuildmodel(m_BFieldMod, m_BFieldModRep);

        // alglib::rbfsetalgohierarchical(m_BFieldPhase, m_PointsDist, 7, 0.0);
        alglib::rbfsetalgohierarchical(m_BFieldPhase, m_PointsDist/4, 3, 0.0);
        // alglib::rbfsetalgobiharmonic(m_BFieldPhase, 0.0);
        // alglib::rbfsetalgomultiquadricauto(m_BFieldPhase, 0.0);
        // alglib::rbfsetalgothinplatespline(m_BFieldPhase, 0.0);
        alglib::rbfbuildmodel(m_BFieldPhase, m_BFieldPhaseRep);
        

        std::cout << "EMod check (should be 1) ---> " << int(m_EFieldModRep.terminationtype) << std::endl;
        std::cout << "EPhase check (should be 1) ---> " << int(m_EFieldPhaseRep.terminationtype) << std::endl;
        std::cout << "BMod check (should be 1) ---> " << int(m_BFieldModRep.terminationtype) << std::endl;
        std::cout << "BPhase check (should be 1) ---> " << int(m_BFieldPhaseRep.terminationtype) << std::endl;
    }

    FieldValue get(const Position& pos) const override
    {
        FieldValue result;
        // alglib::real_1d_array coord, EReal, EImag, BReal, BImag;
        alglib::real_1d_array coord, EMod, EPhase, BMod, BPhase;
        coord.setlength(3);
        // EReal.setlength(3);
        // EImag.setlength(3);
        // BReal.setlength(3);
        // BImag.setlength(3);
        EMod.setlength(3);
        EPhase.setlength(3);
        BMod.setlength(3);
        BPhase.setlength(3);
        coord[0]=pos[0];
        coord[1]=pos[1];
        coord[2]=pos[2];

        alglib::rbfcalcbuffer ERBUF, EIBUF, BRBUF, BIBUF;

        // alglib::rbfcreatecalcbuffer(m_EFieldReal, ERBUF);
        // alglib::rbfcreatecalcbuffer(m_EFieldImag, EIBUF);
        // alglib::rbfcreatecalcbuffer(m_BFieldReal, BRBUF);
        // alglib::rbfcreatecalcbuffer(m_BFieldImag, BIBUF);

        alglib::rbfcreatecalcbuffer(m_EFieldMod, ERBUF);
        alglib::rbfcreatecalcbuffer(m_EFieldPhase, EIBUF);
        alglib::rbfcreatecalcbuffer(m_BFieldMod, BRBUF);
        alglib::rbfcreatecalcbuffer(m_BFieldPhase, BIBUF);
    
        // alglib::rbftscalcbuf(m_EFieldReal, ERBUF, coord, EReal);
        // alglib::rbftscalcbuf(m_EFieldImag, EIBUF, coord, EImag);
        // alglib::rbftscalcbuf(m_BFieldReal, BRBUF, coord, BReal);
        // alglib::rbftscalcbuf(m_BFieldImag, BIBUF, coord, BImag);

        alglib::rbftscalcbuf(m_EFieldMod, ERBUF, coord, EMod);
        alglib::rbftscalcbuf(m_EFieldPhase, EIBUF, coord, EPhase);
        alglib::rbftscalcbuf(m_BFieldMod, BRBUF, coord, BMod);
        alglib::rbftscalcbuf(m_BFieldPhase, BIBUF, coord, BPhase);



        // result.E[0] = EReal[0]+std::complex<double>(0.0, 1.0)*EImag[0];
        // result.E[1] = EReal[1]+std::complex<double>(0.0, 1.0)*EImag[1];
        // result.E[2] = EReal[2]+std::complex<double>(0.0, 1.0)*EImag[2];
        // result.B[0] = BReal[0]+std::complex<double>(0.0, 1.0)*BImag[0];
        // result.B[1] = BReal[1]+std::complex<double>(0.0, 1.0)*BImag[1];
        // result.B[2] = BReal[2]+std::complex<double>(0.0, 1.0)*BImag[2];

        result.E[0] = EMod[0]*exp(std::complex<double>(0.0, 1.0)*EPhase[0]);
        result.E[1] = EMod[1]*exp(std::complex<double>(0.0, 1.0)*EPhase[1]);
        result.E[2] = EMod[2]*exp(std::complex<double>(0.0, 1.0)*EPhase[2]);
        result.B[0] = BMod[0]*exp(std::complex<double>(0.0, 1.0)*BPhase[0]);
        result.B[1] = BMod[1]*exp(std::complex<double>(0.0, 1.0)*BPhase[1]);
        result.B[2] = BMod[2]*exp(std::complex<double>(0.0, 1.0)*BPhase[2]);

        return result;
    }

private:
    size_t m_npoints;
    double m_PointsDist;
    std::array<double, n_points> m_xs, m_ys, m_zs;
    std::array<FieldValue, n_points> &m_FieldInRegion;
    // alglib::rbfmodel m_EFieldReal, m_EFieldImag, m_BFieldReal, m_BFieldImag;
    // alglib::rbfreport m_EFieldRealRep, m_EFieldImagRep, m_BFieldRealRep, m_BFieldImagRep;
    // alglib::real_2d_array m_Points_and_EReals, m_Points_and_EImags, m_Points_and_BReals, m_Points_and_BImags;
    alglib::rbfmodel m_EFieldMod, m_EFieldPhase, m_BFieldMod, m_BFieldPhase;
    alglib::rbfreport m_EFieldModRep, m_EFieldPhaseRep, m_BFieldModRep, m_BFieldPhaseRep;
    alglib::real_2d_array m_Points_and_EMods, m_Points_and_EPhases, m_Points_and_BMods, m_Points_and_BPhases;
    // std::array<std::array<double, 6>, n_points> ;
};

/*
template <const size_t n_points>
class FieldInterpolated1 : public FieldBase
{
public:
    FieldInterpolated1(double lambda, std::array<double, n_points> &xs, std::array<double, n_points> &ys, std::array<double, n_points> &zs, double PointsDist, std::array<FieldValue, n_points> &FieldInRegion) :
        FieldBase(lambda), m_npoints(n_points), m_xs(xs),  m_ys(ys),  m_zs(zs), m_FieldInRegion(FieldInRegion), m_PointsDist(PointsDist)
    {

        alglib::idwbuildercreate(3,3, m_EFieldModBuild);
        alglib::idwbuildercreate(3,3, m_EFieldPhaseBuild);
        alglib::idwbuildercreate(3,3, m_BFieldModBuild);
        alglib::idwbuildercreate(3,3, m_BFieldPhaseBuild);

        m_Points_and_EMods.setlength(m_npoints, 6);
        m_Points_and_EPhases.setlength(m_npoints, 6);
        m_Points_and_BMods.setlength(m_npoints, 6);
        m_Points_and_BPhases.setlength(m_npoints, 6);
       
        for (int i=0; i<m_npoints; i++) {
            m_Points_and_EMods[i][0]=m_xs[i];
            m_Points_and_EMods[i][1]=m_ys[i];
            m_Points_and_EMods[i][2]=m_zs[i];
            m_Points_and_EMods[i][3]=vec_modulus(m_FieldInRegion[i].E)[0];
            m_Points_and_EMods[i][4]=vec_modulus(m_FieldInRegion[i].E)[1];
            m_Points_and_EMods[i][5]=vec_modulus(m_FieldInRegion[i].E)[2];

            m_Points_and_EPhases[i][0]=m_xs[i];
            m_Points_and_EPhases[i][1]=m_ys[i];
            m_Points_and_EPhases[i][2]=m_zs[i];
            m_Points_and_EPhases[i][3]=vec_phases(m_FieldInRegion[i].E)[0];
            m_Points_and_EPhases[i][4]=vec_phases(m_FieldInRegion[i].E)[1];
            m_Points_and_EPhases[i][5]=vec_phases(m_FieldInRegion[i].E)[2];

            m_Points_and_BMods[i][0]=m_xs[i];
            m_Points_and_BMods[i][1]=m_ys[i];
            m_Points_and_BMods[i][2]=m_zs[i];
            m_Points_and_BMods[i][3]=vec_modulus(m_FieldInRegion[i].B)[0];
            m_Points_and_BMods[i][4]=vec_modulus(m_FieldInRegion[i].B)[1];
            m_Points_and_BMods[i][5]=vec_modulus(m_FieldInRegion[i].B)[2];

            m_Points_and_BPhases[i][0]=m_xs[i];
            m_Points_and_BPhases[i][1]=m_ys[i];
            m_Points_and_BPhases[i][2]=m_zs[i];
            m_Points_and_BPhases[i][3]=vec_phases(m_FieldInRegion[i].B)[0];
            m_Points_and_BPhases[i][4]=vec_phases(m_FieldInRegion[i].B)[1];
            m_Points_and_BPhases[i][5]=vec_phases(m_FieldInRegion[i].B)[2];
        }

        alglib::idwbuildersetpoints(m_EFieldModBuild, m_Points_and_EMods);
        alglib::idwbuildersetpoints(m_EFieldPhaseBuild, m_Points_and_EPhases);
        alglib::idwbuildersetpoints(m_BFieldModBuild, m_Points_and_BMods);
        alglib::idwbuildersetpoints(m_BFieldPhaseBuild, m_Points_and_BPhases);


        alglib::idwbuildersetalgomstab(m_EFieldModBuild, m_PointsDist*16);
        alglib::idwfit(m_EFieldModBuild, m_EFieldMod, m_EFieldModRep);

        alglib::idwbuildersetalgomstab(m_EFieldPhaseBuild, m_PointsDist*16);
        alglib::idwfit(m_EFieldPhaseBuild, m_EFieldPhase, m_EFieldPhaseRep);

        alglib::idwbuildersetalgomstab(m_BFieldModBuild, m_PointsDist*16);
        alglib::idwfit(m_BFieldModBuild, m_BFieldMod, m_BFieldModRep);

        alglib::idwbuildersetalgomstab(m_BFieldPhaseBuild, m_PointsDist*16);
        alglib::idwfit(m_BFieldPhaseBuild, m_BFieldPhase, m_BFieldPhaseRep);


        // std::cout << "EMod check (should be 1) ---> " << int(m_EFieldModRep.terminationtype) << std::endl;
        // std::cout << "EPhase check (should be 1) ---> " << int(m_EFieldPhaseRep.terminationtype) << std::endl;
        // std::cout << "BMod check (should be 1) ---> " << int(m_BFieldModRep.terminationtype) << std::endl;
        // std::cout << "BPhase check (should be 1) ---> " << int(m_BFieldPhaseRep.terminationtype) << std::endl;
    }

    FieldValue get(const Position& pos) const override
    {
        FieldValue result;
        alglib::real_1d_array coord, EMod, EPhase, BMod, BPhase;
        coord.setlength(3);
        EMod.setlength(3);
        EPhase.setlength(3);
        BMod.setlength(3);
        BPhase.setlength(3);
        coord[0]=pos[0];
        coord[1]=pos[1];
        coord[2]=pos[2];

        alglib::idwcalcbuffer ERBUF, EIBUF, BRBUF, BIBUF;

        alglib::idwcreatecalcbuffer(m_EFieldMod, ERBUF);
        alglib::idwcreatecalcbuffer(m_EFieldPhase, EIBUF);
        alglib::idwcreatecalcbuffer(m_BFieldMod, BRBUF);
        alglib::idwcreatecalcbuffer(m_BFieldPhase, BIBUF);
    
        alglib::idwtscalcbuf(m_EFieldMod, ERBUF, coord, EMod);
        alglib::idwtscalcbuf(m_EFieldPhase, EIBUF, coord, EPhase);
        alglib::idwtscalcbuf(m_BFieldMod, BRBUF, coord, BMod);
        alglib::idwtscalcbuf(m_BFieldPhase, BIBUF, coord, BPhase);


        result.E[0] = EMod[0]*exp(std::complex<double>(0.0, 1.0)*EPhase[0]);
        result.E[1] = EMod[1]*exp(std::complex<double>(0.0, 1.0)*EPhase[1]);
        result.E[2] = EMod[2]*exp(std::complex<double>(0.0, 1.0)*EPhase[2]);
        result.B[0] = BMod[0]*exp(std::complex<double>(0.0, 1.0)*BPhase[0]);
        result.B[1] = BMod[1]*exp(std::complex<double>(0.0, 1.0)*BPhase[1]);
        result.B[2] = BMod[2]*exp(std::complex<double>(0.0, 1.0)*BPhase[2]);

        return result;
    }

private:
    size_t m_npoints;
    double m_PointsDist;
    std::array<double, n_points> m_xs, m_ys, m_zs;
    std::array<FieldValue, n_points> &m_FieldInRegion;
    // alglib::rbfmodel m_EFieldReal, m_EFieldImag, m_BFieldReal, m_BFieldImag;
    // alglib::rbfreport m_EFieldRealRep, m_EFieldImagRep, m_BFieldRealRep, m_BFieldImagRep;
    // alglib::real_2d_array m_Points_and_EReals, m_Points_and_EImags, m_Points_and_BReals, m_Points_and_BImags;
    alglib::idwbuilder m_EFieldModBuild, m_EFieldPhaseBuild, m_BFieldModBuild, m_BFieldPhaseBuild;
    alglib::idwmodel m_EFieldMod, m_EFieldPhase, m_BFieldMod, m_BFieldPhase;
    alglib::idwreport m_EFieldModRep, m_EFieldPhaseRep, m_BFieldModRep, m_BFieldPhaseRep;
    alglib::real_2d_array m_Points_and_EMods, m_Points_and_EPhases, m_Points_and_BMods, m_Points_and_BPhases;
    // std::array<std::array<double, 6>, n_points> ;
};
*/


template <const size_t n_points>
class FieldInterpolated2d : public FieldBase
{
public:
    FieldInterpolated2d(double lambda, std::array<double, n_points> &xs, std::array<double, n_points> &ys, double PointsDist, std::array<std::array<FieldValue, n_points>, n_points> &FieldInRegion, ISurface& surf) :
        FieldBase(lambda), m_npoints(n_points), m_xs(xs),  m_ys(ys),  m_FieldInRegion(FieldInRegion), m_PointsDist(PointsDist), m_surf(surf)
    {

        // alglib::rbfcreate(2,3, m_EFieldReal);
        // alglib::rbfcreate(2,3, m_EFieldImag);
        // alglib::rbfcreate(2,3, m_BFieldReal);
        // alglib::rbfcreate(2,3, m_BFieldImag);

        // m_Points_and_EReals.setlength(m_npoints*m_npoints, 5);
        // m_Points_and_EImags.setlength(m_npoints*m_npoints, 5);
        // m_Points_and_BReals.setlength(m_npoints*m_npoints, 5);
        // m_Points_and_BImags.setlength(m_npoints*m_npoints, 5);
     
        // size_t count;
       
        // for (int i=0; i<m_npoints; i++) {
        //     for (int j=0; j<m_npoints; j++) {

        //         count = i*m_npoints + j;

        //         m_Points_and_EReals[count][0]=m_xs[i];
        //         m_Points_and_EReals[count][1]=m_ys[j];
        //         m_Points_and_EReals[count][2]=vec_real(m_FieldInRegion[i][j].E)[0];
        //         m_Points_and_EReals[count][3]=vec_real(m_FieldInRegion[i][j].E)[1];
        //         m_Points_and_EReals[count][4]=vec_real(m_FieldInRegion[i][j].E)[2];

        //         m_Points_and_EImags[count][0]=m_xs[i];
        //         m_Points_and_EImags[count][1]=m_ys[j];
        //         m_Points_and_EImags[count][2]=vec_imag(m_FieldInRegion[i][j].E)[0];
        //         m_Points_and_EImags[count][3]=vec_imag(m_FieldInRegion[i][j].E)[1];
        //         m_Points_and_EImags[count][4]=vec_imag(m_FieldInRegion[i][j].E)[2];

        //         m_Points_and_BReals[count][0]=m_xs[i];
        //         m_Points_and_BReals[count][1]=m_ys[j];
        //         m_Points_and_BReals[count][2]=vec_real(m_FieldInRegion[i][j].B)[0];
        //         m_Points_and_BReals[count][3]=vec_real(m_FieldInRegion[i][j].B)[1];
        //         m_Points_and_BReals[count][4]=vec_real(m_FieldInRegion[i][j].B)[2];

        //         m_Points_and_BImags[count][0]=m_xs[i];
        //         m_Points_and_BImags[count][1]=m_ys[j];
        //         m_Points_and_BImags[count][2]=vec_imag(m_FieldInRegion[i][j].B)[0];
        //         m_Points_and_BImags[count][3]=vec_imag(m_FieldInRegion[i][j].B)[1];
        //         m_Points_and_BImags[count][4]=vec_imag(m_FieldInRegion[i][j].B)[2];
        //     }
        // }

        // alglib::rbfsetpoints(m_EFieldReal, m_Points_and_EReals);
        // alglib::rbfsetpoints(m_EFieldImag, m_Points_and_EImags);
        // alglib::rbfsetpoints(m_BFieldReal, m_Points_and_BReals);
        // alglib::rbfsetpoints(m_BFieldImag, m_Points_and_BImags);

        alglib::rbfcreate(2,3, m_EFieldMod);
        alglib::rbfcreate(2,3, m_EFieldPhase);
        alglib::rbfcreate(2,3, m_BFieldMod);
        alglib::rbfcreate(2,3, m_BFieldPhase);

        m_Points_and_EMods.setlength(m_npoints*m_npoints, 5);
        m_Points_and_EPhases.setlength(m_npoints*m_npoints, 5);
        m_Points_and_BMods.setlength(m_npoints*m_npoints, 5);
        m_Points_and_BPhases.setlength(m_npoints*m_npoints, 5);

        size_t count;
       
        for (int i=0; i<m_npoints; i++) {
            for (int j=0; j<m_npoints; j++) {

                count = i*m_npoints + j;

                m_Points_and_EMods[count][0]=m_xs[i];
                m_Points_and_EMods[count][1]=m_ys[j];
                m_Points_and_EMods[count][2]=vec_modulus(m_FieldInRegion[i][j].E)[0];
                m_Points_and_EMods[count][3]=vec_modulus(m_FieldInRegion[i][j].E)[1];
                m_Points_and_EMods[count][4]=vec_modulus(m_FieldInRegion[i][j].E)[2];

                m_Points_and_EPhases[count][0]=m_xs[i];
                m_Points_and_EPhases[count][1]=m_ys[j];
                m_Points_and_EPhases[count][2]=vec_phases(m_FieldInRegion[i][j].E)[0];
                m_Points_and_EPhases[count][3]=vec_phases(m_FieldInRegion[i][j].E)[1];
                m_Points_and_EPhases[count][4]=vec_phases(m_FieldInRegion[i][j].E)[2];

                m_Points_and_BMods[count][0]=m_xs[i];
                m_Points_and_BMods[count][1]=m_ys[j];
                m_Points_and_BMods[count][2]=vec_modulus(m_FieldInRegion[i][j].B)[0];
                m_Points_and_BMods[count][3]=vec_modulus(m_FieldInRegion[i][j].B)[1];
                m_Points_and_BMods[count][4]=vec_modulus(m_FieldInRegion[i][j].B)[2];

                m_Points_and_BPhases[count][0]=m_xs[i];
                m_Points_and_BPhases[count][1]=m_ys[j];
                m_Points_and_BPhases[count][2]=vec_phases(m_FieldInRegion[i][j].B)[0];
                m_Points_and_BPhases[count][3]=vec_phases(m_FieldInRegion[i][j].B)[1];
                m_Points_and_BPhases[count][4]=vec_phases(m_FieldInRegion[i][j].B)[2];
            }
        }

        alglib::rbfsetpoints(m_EFieldMod, m_Points_and_EMods);
        alglib::rbfsetpoints(m_EFieldPhase, m_Points_and_EPhases);
        alglib::rbfsetpoints(m_BFieldMod, m_Points_and_BMods);
        alglib::rbfsetpoints(m_BFieldPhase, m_Points_and_BPhases);

        struct Task{
            alglib::rbfmodel &m_EFieldMod, &m_EFieldPhase, &m_BFieldMod, &m_BFieldPhase;
            alglib::rbfreport &m_EFieldModRep, &m_EFieldPhaseRep, &m_BFieldModRep, &m_BFieldPhaseRep;
            double m_PointsDist;
            size_t j;
            void operator()(){
                if (j==0) {
                    alglib::rbfsetalgohierarchical(m_EFieldMod, m_PointsDist, 7, 1.0E-4);
                    alglib::rbfbuildmodel(m_EFieldMod, m_EFieldModRep);
                    std::cout << "EMod check (should be 1) ---> " << int(m_EFieldModRep.terminationtype) << std::endl;
                };

                if (j==1) {
                    alglib::rbfsetalgohierarchical(m_EFieldPhase, m_PointsDist, 7, 1.0E-4);
                    alglib::rbfbuildmodel(m_EFieldPhase, m_EFieldPhaseRep);
                    std::cout << "EPhase check (should be 1) ---> " << int(m_EFieldPhaseRep.terminationtype) << std::endl;
                };

                if (j==2) {
                    alglib::rbfsetalgohierarchical(m_BFieldMod, m_PointsDist, 7, 1.0E-4);
                    alglib::rbfbuildmodel(m_BFieldMod, m_BFieldModRep);
                    std::cout << "BMod check (should be 1) ---> " << int(m_BFieldModRep.terminationtype) << std::endl;
                };

                if (j==3) {
                    alglib::rbfsetalgohierarchical(m_BFieldPhase, m_PointsDist, 7, 1.0E-4);
                    alglib::rbfbuildmodel(m_BFieldPhase, m_BFieldPhaseRep);
                    std::cout << "BPhase check (should be 1) ---> " << int(m_BFieldPhaseRep.terminationtype) << std::endl;
                };
                
            }
            operator bool()const{
                return true;
            }
        };
        using Threads = Threading::ThreadPool<ConcurrentQueue<Task>,4>;
        typename Threads::TaskQueue tasks;

        for (size_t j = 0; j < 4; j++)
        {
            tasks.push(Task{m_EFieldMod, m_EFieldPhase, m_BFieldMod, m_BFieldPhase, m_EFieldModRep, m_EFieldPhaseRep, m_BFieldModRep, m_BFieldPhaseRep, m_PointsDist, j});
        }
        Threads threads{tasks};
        threads.run();

        // struct Task{
        //     alglib::rbfmodel m_EFieldReal, m_EFieldImag, m_BFieldReal, m_BFieldImag;
        //     alglib::rbfreport m_EFieldRealRep, m_EFieldImagRep, m_BFieldRealRep, m_BFieldImagRep;
        //     double m_PointsDist;
        //     size_t j;
        //     void operator()(){
        //         if (j==0) {
        //             alglib::rbfsetalgohierarchical(m_EFieldReal, m_PointsDist, 7, 1.0E-4);
        //             alglib::rbfbuildmodel(m_EFieldReal, m_EFieldRealRep);
        //             std::cout << "EReal check (should be 1) ---> " << int(m_EFieldRealRep.terminationtype) << std::endl;
        //         };

        //         if (j==1) {
        //             alglib::rbfsetalgohierarchical(m_EFieldImag, m_PointsDist, 7, 1.0E-4);
        //             alglib::rbfbuildmodel(m_EFieldImag, m_EFieldImagRep);
        //             std::cout << "EImag check (should be 1) ---> " << int(m_EFieldImagRep.terminationtype) << std::endl;
        //         };

        //         if (j==2) {
        //             alglib::rbfsetalgohierarchical(m_BFieldReal, m_PointsDist, 7, 1.0E-4);
        //             alglib::rbfbuildmodel(m_BFieldReal, m_BFieldRealRep);
        //             std::cout << "BReal check (should be 1) ---> " << int(m_BFieldRealRep.terminationtype) << std::endl;
        //         };

        //         if (j==3) {
        //             alglib::rbfsetalgohierarchical(m_BFieldImag, m_PointsDist, 7, 1.0E-4);
        //             alglib::rbfbuildmodel(m_BFieldImag, m_BFieldImagRep);
        //             std::cout << "BImag check (should be 1) ---> " << int(m_BFieldImagRep.terminationtype) << std::endl;
        //         };
                
        //     }
        //     operator bool()const{
        //         return true;
        //     }
        // };
        // using Threads = Threading::ThreadPool<ConcurrentQueue<Task>,4>;
        // typename Threads::TaskQueue tasks;

        // for (size_t j = 0; j < 4; j++)
        // {
        //     tasks.push(Task{m_EFieldReal, m_EFieldImag, m_BFieldReal, m_BFieldImag, m_EFieldRealRep, m_EFieldImagRep, m_BFieldRealRep, m_BFieldImagRep, m_PointsDist, j});
        // }
        // Threads threads{tasks};
        // threads.run();

        

        // alglib::rbfsetalgohierarchical(m_EFieldReal, m_PointsDist, 6, 1.0E-6);
        // // alglib::rbfsetalgobiharmonic(m_EFieldReal, 1.0E-3);
        // // alglib::rbfsetalgomultiquadricauto(m_EFieldReal, 1.0E-3);
        // alglib::rbfbuildmodel(m_EFieldReal, m_EFieldRealRep);

        // alglib::rbfsetalgohierarchical(m_EFieldImag, m_PointsDist, 6, 1.0E-6);
        // // alglib::rbfsetalgobiharmonic(m_EFieldImag, 1.0E-3);
        // // alglib::rbfsetalgomultiquadricauto(m_EFieldImag, 1.0E-3);
        // alglib::rbfbuildmodel(m_EFieldImag, m_EFieldImagRep);

        // alglib::rbfsetalgohierarchical(m_BFieldReal, m_PointsDist, 6, 1.0E-6);
        // // alglib::rbfsetalgobiharmonic(m_BFieldReal, 1.0E-3);
        // // alglib::rbfsetalgomultiquadricauto(m_BFieldReal, 1.0E-3);
        // alglib::rbfbuildmodel(m_BFieldReal, m_BFieldRealRep);

        // alglib::rbfsetalgohierarchical(m_BFieldImag, m_PointsDist, 6, 1.0E-6);
        // // alglib::rbfsetalgobiharmonic(m_BFieldImag, 1.0E-3);
        // // alglib::rbfsetalgomultiquadricauto(m_BFieldImag, 1.0E-3);
        // alglib::rbfbuildmodel(m_BFieldImag, m_BFieldImagRep);
        

        // std::cout << "EReal check (should be 1) ---> " << int(m_EFieldRealRep.terminationtype) << std::endl;
        // std::cout << "EImag check (should be 1) ---> " << int(m_EFieldImagRep.terminationtype) << std::endl;
        // std::cout << "BReal check (should be 1) ---> " << int(m_BFieldRealRep.terminationtype) << std::endl;
        // std::cout << "BImag check (should be 1) ---> " << int(m_BFieldImagRep.terminationtype) << std::endl;


        // alglib::rbfsetalgohierarchical(m_EFieldMod, m_PointsDist, 7, 0.0);
        // // alglib::rbfsetalgobiharmonic(m_EFieldMod, 0.0);
        // // alglib::rbfsetalgomultiquadricauto(m_EFieldMod, 0.0);
        // // alglib::rbfsetalgothinplatespline(m_EFieldMod, 1.0E-3);
        // alglib::rbfbuildmodel(m_EFieldMod, m_EFieldModRep);

        // alglib::rbfsetalgohierarchical(m_EFieldPhase, m_PointsDist, 7, 0.0);
        // // alglib::rbfsetalgobiharmonic(m_EFieldPhase, 0.0);
        // // alglib::rbfsetalgomultiquadricauto(m_EFieldPhase,0.0);
        // // alglib::rbfsetalgothinplatespline(m_EFieldPhase,0.0);
        // alglib::rbfbuildmodel(m_EFieldPhase, m_EFieldPhaseRep);

        // alglib::rbfsetalgohierarchical(m_BFieldMod, m_PointsDist, 7, 0.0);
        // // alglib::rbfsetalgobiharmonic(m_BFieldMod, 0.0);
        // // alglib::rbfsetalgomultiquadricauto(m_BFieldMod, 0.0);
        // // alglib::rbfsetalgothinplatespline(m_BFieldMod, 0.0);
        // alglib::rbfbuildmodel(m_BFieldMod, m_BFieldModRep);

        // alglib::rbfsetalgohierarchical(m_BFieldPhase, m_PointsDist, 7, 0.0);
        // // alglib::rbfsetalgobiharmonic(m_BFieldPhase, 0.0);
        // // alglib::rbfsetalgomultiquadricauto(m_BFieldPhase, 0.0);
        // // alglib::rbfsetalgothinplatespline(m_BFieldPhase, 0.0);
        // alglib::rbfbuildmodel(m_BFieldPhase, m_BFieldPhaseRep);
        

        // std::cout << "EMod check (should be 1) ---> " << int(m_EFieldModRep.terminationtype) << std::endl;
        // std::cout << "EPhase check (should be 1) ---> " << int(m_EFieldPhaseRep.terminationtype) << std::endl;
        // std::cout << "BMod check (should be 1) ---> " << int(m_BFieldModRep.terminationtype) << std::endl;
        // std::cout << "BPhase check (should be 1) ---> " << int(m_BFieldPhaseRep.terminationtype) << std::endl;
    }

    
    FieldValue get(const Position& pos) const override
    {
        FieldValue result;
        // alglib::real_1d_array coord, EReal, EImag, BReal, BImag;
        alglib::real_1d_array coord, EMod, EPhase, BMod, BPhase;
        coord.setlength(2);
        EMod.setlength(3);
        EPhase.setlength(3);
        BMod.setlength(3);
        BPhase.setlength(3);

        Vector2D pos2d = m_surf.point2d(pos);

        coord[0]=pos2d[0];
        coord[1]=pos2d[1];

        alglib::rbfcalcbuffer ERBUF, EIBUF, BRBUF, BIBUF;

        alglib::rbfcreatecalcbuffer(m_EFieldMod, ERBUF);
        alglib::rbfcreatecalcbuffer(m_EFieldPhase, EIBUF);
        alglib::rbfcreatecalcbuffer(m_BFieldMod, BRBUF);
        alglib::rbfcreatecalcbuffer(m_BFieldPhase, BIBUF);

        alglib::rbftscalcbuf(m_EFieldMod, ERBUF, coord, EMod);
        alglib::rbftscalcbuf(m_EFieldPhase, EIBUF, coord, EPhase);
        alglib::rbftscalcbuf(m_BFieldMod, BRBUF, coord, BMod);
        alglib::rbftscalcbuf(m_BFieldPhase, BIBUF, coord, BPhase);

        result.E[0] = EMod[0]*exp(std::complex<double>(0.0, 1.0)*EPhase[0]);
        result.E[1] = EMod[1]*exp(std::complex<double>(0.0, 1.0)*EPhase[1]);
        result.E[2] = EMod[2]*exp(std::complex<double>(0.0, 1.0)*EPhase[2]);
        result.B[0] = BMod[0]*exp(std::complex<double>(0.0, 1.0)*BPhase[0]);
        result.B[1] = BMod[1]*exp(std::complex<double>(0.0, 1.0)*BPhase[1]);
        result.B[2] = BMod[2]*exp(std::complex<double>(0.0, 1.0)*BPhase[2]);

        // EReal.setlength(3);
        // EImag.setlength(3);
        // BReal.setlength(3);
        // BImag.setlength(3);

        // alglib::rbfcreatecalcbuffer(m_EFieldReal, ERBUF);
        // alglib::rbfcreatecalcbuffer(m_EFieldImag, EIBUF);
        // alglib::rbfcreatecalcbuffer(m_BFieldReal, BRBUF);
        // alglib::rbfcreatecalcbuffer(m_BFieldImag, BIBUF);
    
        // alglib::rbftscalcbuf(m_EFieldReal, ERBUF, coord, EReal);
        // alglib::rbftscalcbuf(m_EFieldImag, EIBUF, coord, EImag);
        // alglib::rbftscalcbuf(m_BFieldReal, BRBUF, coord, BReal);
        // alglib::rbftscalcbuf(m_BFieldImag, BIBUF, coord, BImag);

        // result.E[0] = EReal[0]+std::complex<double>(0.0, 1.0)*EImag[0];
        // result.E[1] = EReal[1]+std::complex<double>(0.0, 1.0)*EImag[1];
        // result.E[2] = EReal[2]+std::complex<double>(0.0, 1.0)*EImag[2];
        // result.B[0] = BReal[0]+std::complex<double>(0.0, 1.0)*BImag[0];
        // result.B[1] = BReal[1]+std::complex<double>(0.0, 1.0)*BImag[1];
        // result.B[2] = BReal[2]+std::complex<double>(0.0, 1.0)*BImag[2];

        return result;
    }

private:
    size_t m_npoints;
    double m_PointsDist;
    ISurface& m_surf;
    std::array<double, n_points> m_xs, m_ys;
    std::array<std::array<FieldValue, n_points>, n_points> &m_FieldInRegion;
    // alglib::rbfmodel m_EFieldReal, m_EFieldImag, m_BFieldReal, m_BFieldImag;
    // alglib::rbfreport m_EFieldRealRep, m_EFieldImagRep, m_BFieldRealRep, m_BFieldImagRep;
    // alglib::real_2d_array m_Points_and_EReals, m_Points_and_EImags, m_Points_and_BReals, m_Points_and_BImags;
    alglib::rbfmodel m_EFieldMod, m_EFieldPhase, m_BFieldMod, m_BFieldPhase;
    alglib::rbfreport m_EFieldModRep, m_EFieldPhaseRep, m_BFieldModRep, m_BFieldPhaseRep;
    alglib::real_2d_array m_Points_and_EMods, m_Points_and_EPhases, m_Points_and_BMods, m_Points_and_BPhases;
    // std::array<std::array<double, 6>, n_points> ;
};


/*
template <const size_t n_points>
class FieldInterpolatedMPComps : public FieldBase
{
public:
    FieldInterpolatedMPComps(double lambda, std::array<double, n_points> &xs, std::array<double, n_points> &ys, std::array<double, n_points> &zs, double PointsDist, std::array<FieldValue, n_points> &FieldInRegion) :
        FieldBase(lambda), m_npoints(n_points), m_xs(xs),  m_ys(ys),  m_zs(zs), m_FieldInRegion(FieldInRegion), m_PointsDist(PointsDist)

    {
        alglib::rbfcreate(3,1, m_EFieldModx);
        alglib::rbfcreate(3,1, m_EFieldMody);
        alglib::rbfcreate(3,1, m_EFieldModz);
        alglib::rbfcreate(3,1, m_BFieldModx);
        alglib::rbfcreate(3,1, m_BFieldMody);
        alglib::rbfcreate(3,1, m_BFieldModz);
        alglib::rbfcreate(3,1, m_EFieldPhasex);
        alglib::rbfcreate(3,1, m_EFieldPhasey);
        alglib::rbfcreate(3,1, m_EFieldPhasez);
        alglib::rbfcreate(3,1, m_BFieldPhasex);
        alglib::rbfcreate(3,1, m_BFieldPhasey);
        alglib::rbfcreate(3,1, m_BFieldPhasez);

        m_Points_and_EModx.setlength(m_npoints, 4);
        m_Points_and_EMody.setlength(m_npoints, 4);
        m_Points_and_EModz.setlength(m_npoints, 4);
        m_Points_and_EPhasex.setlength(m_npoints, 4);
        m_Points_and_EPhasey.setlength(m_npoints, 4);
        m_Points_and_EPhasez.setlength(m_npoints, 4);
        m_Points_and_BModx.setlength(m_npoints, 4);
        m_Points_and_BMody.setlength(m_npoints, 4);
        m_Points_and_BModz.setlength(m_npoints, 4);
        m_Points_and_BPhasex.setlength(m_npoints, 4);
        m_Points_and_BPhasey.setlength(m_npoints, 4);
        m_Points_and_BPhasez.setlength(m_npoints, 4);

       
        for (int i=0; i<m_npoints; i++) {
            m_Points_and_EModx[i][0]=m_xs[i];
            m_Points_and_EModx[i][1]=m_ys[i];
            m_Points_and_EModx[i][2]=m_zs[i];
            m_Points_and_EModx[i][3]=vec_modulus(m_FieldInRegion[i].E)[0];

            m_Points_and_EMody[i][0]=m_xs[i];
            m_Points_and_EMody[i][1]=m_ys[i];
            m_Points_and_EMody[i][2]=m_zs[i];
            m_Points_and_EMody[i][3]=vec_modulus(m_FieldInRegion[i].E)[1];

            m_Points_and_EModz[i][0]=m_xs[i];
            m_Points_and_EModz[i][1]=m_ys[i];
            m_Points_and_EModz[i][2]=m_zs[i];
            m_Points_and_EModz[i][3]=vec_modulus(m_FieldInRegion[i].E)[2];

            m_Points_and_EPhasex[i][0]=m_xs[i];
            m_Points_and_EPhasex[i][1]=m_ys[i];
            m_Points_and_EPhasex[i][2]=m_zs[i];
            m_Points_and_EPhasex[i][3]=vec_phases(m_FieldInRegion[i].E)[0];

            m_Points_and_EPhasey[i][0]=m_xs[i];
            m_Points_and_EPhasey[i][1]=m_ys[i];
            m_Points_and_EPhasey[i][2]=m_zs[i];
            m_Points_and_EPhasey[i][3]=vec_phases(m_FieldInRegion[i].E)[1];

            m_Points_and_EPhasez[i][0]=m_xs[i];
            m_Points_and_EPhasez[i][1]=m_ys[i];
            m_Points_and_EPhasez[i][2]=m_zs[i];
            m_Points_and_EPhasez[i][3]=vec_phases(m_FieldInRegion[i].E)[2];

            m_Points_and_BModx[i][0]=m_xs[i];
            m_Points_and_BModx[i][1]=m_ys[i];
            m_Points_and_BModx[i][2]=m_zs[i];
            m_Points_and_BModx[i][3]=vec_modulus(m_FieldInRegion[i].B)[0];

            m_Points_and_BMody[i][0]=m_xs[i];
            m_Points_and_BMody[i][1]=m_ys[i];
            m_Points_and_BMody[i][2]=m_zs[i];
            m_Points_and_BMody[i][3]=vec_modulus(m_FieldInRegion[i].B)[1];

            m_Points_and_BModz[i][0]=m_xs[i];
            m_Points_and_BModz[i][1]=m_ys[i];
            m_Points_and_BModz[i][2]=m_zs[i];
            m_Points_and_BModz[i][3]=vec_modulus(m_FieldInRegion[i].B)[2];

            m_Points_and_BPhasex[i][0]=m_xs[i];
            m_Points_and_BPhasex[i][1]=m_ys[i];
            m_Points_and_BPhasex[i][2]=m_zs[i];
            m_Points_and_BPhasex[i][3]=vec_phases(m_FieldInRegion[i].E)[0];

            m_Points_and_BPhasey[i][0]=m_xs[i];
            m_Points_and_BPhasey[i][1]=m_ys[i];
            m_Points_and_BPhasey[i][2]=m_zs[i];
            m_Points_and_BPhasey[i][3]=vec_phases(m_FieldInRegion[i].E)[1];

            m_Points_and_BPhasez[i][0]=m_xs[i];
            m_Points_and_BPhasez[i][1]=m_ys[i];
            m_Points_and_BPhasez[i][2]=m_zs[i];
            m_Points_and_BPhasez[i][3]=vec_phases(m_FieldInRegion[i].E)[2];
        }

        alglib::rbfsetpoints(m_EFieldModx, m_Points_and_EModx);
        alglib::rbfsetpoints(m_EFieldMody, m_Points_and_EMody);
        alglib::rbfsetpoints(m_EFieldModz, m_Points_and_EModz);
        alglib::rbfsetpoints(m_EFieldPhasex, m_Points_and_EPhasex);
        alglib::rbfsetpoints(m_EFieldPhasey, m_Points_and_EPhasey);
        alglib::rbfsetpoints(m_EFieldPhasez, m_Points_and_EPhasez);
        alglib::rbfsetpoints(m_BFieldModx, m_Points_and_BModx);
        alglib::rbfsetpoints(m_BFieldMody, m_Points_and_BMody);
        alglib::rbfsetpoints(m_BFieldModz, m_Points_and_BModz);
        alglib::rbfsetpoints(m_BFieldPhasex, m_Points_and_BPhasex);
        alglib::rbfsetpoints(m_BFieldPhasey, m_Points_and_BPhasey);
        alglib::rbfsetpoints(m_BFieldPhasez, m_Points_and_BPhasez);


        // alglib::rbfsetalgohierarchical(m_EFieldModx, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_EFieldModx, 0.0);
        alglib::rbfbuildmodel(m_EFieldModx, m_EFieldModRepx);

        // alglib::rbfsetalgohierarchical(m_EFieldMody, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_EFieldMody, 0.0);
        // alglib::rbfbuildmodel(m_EFieldMody, m_EFieldModRepy);

        // alglib::rbfsetalgohierarchical(m_EFieldModz, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_EFieldModz, 0.0);
        alglib::rbfbuildmodel(m_EFieldModz, m_EFieldModRepz);

        // alglib::rbfsetalgohierarchical(m_EFieldPhasex, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_EFieldPhasex,0.0);
        alglib::rbfbuildmodel(m_EFieldPhasex, m_EFieldPhaseRepx);

        // alglib::rbfsetalgohierarchical(m_EFieldPhasey, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_EFieldPhasey,0.0);
        alglib::rbfbuildmodel(m_EFieldPhasey, m_EFieldPhaseRepy);

        // alglib::rbfsetalgohierarchical(m_EFieldPhasez, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_EFieldPhasez,0.0);
        alglib::rbfbuildmodel(m_EFieldPhasez, m_EFieldPhaseRepz);

        // alglib::rbfsetalgohierarchical(m_BFieldModx, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_BFieldModx, 0.0);
        alglib::rbfbuildmodel(m_BFieldModx, m_BFieldModRepx);

        // alglib::rbfsetalgohierarchical(m_BFieldMody, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_BFieldMody, 0.0);
        alglib::rbfbuildmodel(m_BFieldMody, m_BFieldModRepy);

        // alglib::rbfsetalgohierarchical(m_BFieldModz, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_BFieldModz, 0.0);
        alglib::rbfbuildmodel(m_BFieldModz, m_BFieldModRepz);

        // alglib::rbfsetalgohierarchical(m_BFieldPhasex, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_BFieldPhasex, 0.0);
        alglib::rbfbuildmodel(m_BFieldPhasex, m_BFieldPhaseRepx);

        // alglib::rbfsetalgohierarchical(m_BFieldPhasey, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_BFieldPhasey, 0.0);
        alglib::rbfbuildmodel(m_BFieldPhasey, m_BFieldPhaseRepy);

        // alglib::rbfsetalgohierarchical(m_BFieldPhasez, m_PointsDist, 5, 0.0);
        alglib::rbfsetalgothinplatespline(m_BFieldPhasez, 0.0);
        alglib::rbfbuildmodel(m_BFieldPhasez, m_BFieldPhaseRepz);
        

        std::cout << "EModx check (should be 1) ---> " << int(m_EFieldModRepx.terminationtype) << std::endl;
        std::cout << "EMody check (should be 1) ---> " << int(m_EFieldModRepy.terminationtype) << std::endl;
        std::cout << "EModz check (should be 1) ---> " << int(m_EFieldModRepz.terminationtype) << std::endl;
        std::cout << "EPhasex check (should be 1) ---> " << int(m_EFieldPhaseRepx.terminationtype) << std::endl;
        std::cout << "EPhasey check (should be 1) ---> " << int(m_EFieldPhaseRepy.terminationtype) << std::endl;
        std::cout << "EPhasez check (should be 1) ---> " << int(m_EFieldPhaseRepz.terminationtype) << std::endl;
        std::cout << "BModx check (should be 1) ---> " << int(m_BFieldModRepx.terminationtype) << std::endl;
        std::cout << "BMody check (should be 1) ---> " << int(m_BFieldModRepy.terminationtype) << std::endl;
        std::cout << "BModz check (should be 1) ---> " << int(m_BFieldModRepz.terminationtype) << std::endl;
        std::cout << "BPhasex check (should be 1) ---> " << int(m_BFieldPhaseRepx.terminationtype) << std::endl;
        std::cout << "BPhasey check (should be 1) ---> " << int(m_BFieldPhaseRepy.terminationtype) << std::endl;
        std::cout << "BPhasez check (should be 1) ---> " << int(m_BFieldPhaseRepz.terminationtype) << std::endl;
    }

    
    FieldValue get(const Position& pos) const override
    {
        FieldValue result;
        // alglib::real_1d_array coord, EReal, EImag, BReal, BImag;
        alglib::real_1d_array coord, EModx, EMody, EModz, EPhasex, EPhasey, EPhasez, BModx, BMody, BModz, BPhasex, BPhasey, BPhasez;
        coord.setlength(3);
        EModx.setlength(1);
        EMody.setlength(1);
        EModz.setlength(1);
        EPhasex.setlength(1);
        EPhasey.setlength(1);
        EPhasez.setlength(1);
        BModx.setlength(1);
        BMody.setlength(1);
        BModz.setlength(1);
        BPhasex.setlength(1);
        BPhasey.setlength(1);
        BPhasez.setlength(1);

        coord[0]=pos[0];
        coord[1]=pos[1];
        coord[2]=pos[2];

        alglib::rbfcalcbuffer ERBUFx, ERBUFy, ERBUFz, EIBUFx, EIBUFy, EIBUFz, BRBUFx, BRBUFy, BRBUFz, BIBUFx, BIBUFy, BIBUFz;

        alglib::rbfcreatecalcbuffer(m_EFieldModx, ERBUFx);
        alglib::rbfcreatecalcbuffer(m_EFieldMody, ERBUFy);
        alglib::rbfcreatecalcbuffer(m_EFieldModz, ERBUFz);
        alglib::rbfcreatecalcbuffer(m_EFieldPhasex, EIBUFx);
        alglib::rbfcreatecalcbuffer(m_EFieldPhasey, EIBUFy);
        alglib::rbfcreatecalcbuffer(m_EFieldPhasez, EIBUFz);
        alglib::rbfcreatecalcbuffer(m_BFieldModx, BRBUFx);
        alglib::rbfcreatecalcbuffer(m_BFieldMody, BRBUFy);
        alglib::rbfcreatecalcbuffer(m_BFieldModz, BRBUFz);
        alglib::rbfcreatecalcbuffer(m_BFieldPhasex, BIBUFx);
        alglib::rbfcreatecalcbuffer(m_BFieldPhasey, BIBUFy);
        alglib::rbfcreatecalcbuffer(m_BFieldPhasez, BIBUFz);
        
        alglib::rbftscalcbuf(m_EFieldModx, ERBUFx, coord, EModx);
        alglib::rbftscalcbuf(m_EFieldMody, ERBUFy, coord, EMody);
        alglib::rbftscalcbuf(m_EFieldModz, ERBUFz, coord, EModz);
        alglib::rbftscalcbuf(m_EFieldPhasex, EIBUFx, coord, EPhasex);
        alglib::rbftscalcbuf(m_EFieldPhasey, EIBUFy, coord, EPhasey);
        alglib::rbftscalcbuf(m_EFieldPhasez, EIBUFz, coord, EPhasez);
        alglib::rbftscalcbuf(m_BFieldModx, BRBUFx, coord, BModx);
        alglib::rbftscalcbuf(m_BFieldMody, BRBUFy, coord, BMody);
        alglib::rbftscalcbuf(m_BFieldModz, BRBUFz, coord, BModz);
        alglib::rbftscalcbuf(m_BFieldPhasex, BIBUFx, coord, BPhasex);
        alglib::rbftscalcbuf(m_BFieldPhasey, BIBUFy, coord, BPhasey);
        alglib::rbftscalcbuf(m_BFieldPhasez, BIBUFz, coord, BPhasez);

        result.E[0] = EModx[0]*exp(std::complex<double>(0.0, 1.0)*EPhasex[0]);
        result.E[1] = EMody[0]*exp(std::complex<double>(0.0, 1.0)*EPhasey[0]);
        result.E[2] = EModz[0]*exp(std::complex<double>(0.0, 1.0)*EPhasez[0]);
        result.B[0] = BModx[0]*exp(std::complex<double>(0.0, 1.0)*BPhasex[0]);
        result.B[1] = BMody[0]*exp(std::complex<double>(0.0, 1.0)*BPhasey[0]);
        result.B[2] = BModz[0]*exp(std::complex<double>(0.0, 1.0)*BPhasez[0]);

        return result;
    }

private:
    size_t m_npoints;
    double m_PointsDist;
    std::array<double, n_points> m_xs, m_ys, m_zs;
    std::array<FieldValue, n_points> &m_FieldInRegion;
    // alglib::rbfmodel m_EFieldReal, m_EFieldImag, m_BFieldReal, m_BFieldImag;
    // alglib::rbfreport m_EFieldRealRep, m_EFieldImagRep, m_BFieldRealRep, m_BFieldImagRep;
    // alglib::real_2d_array m_Points_and_EReals, m_Points_and_EImags, m_Points_and_BReals, m_Points_and_BImags;
    alglib::rbfmodel m_EFieldModx, m_EFieldMody, m_EFieldModz, m_EFieldPhasex, m_EFieldPhasey, m_EFieldPhasez, m_BFieldModx, m_BFieldMody, m_BFieldModz, m_BFieldPhasex, m_BFieldPhasey, m_BFieldPhasez;
    alglib::rbfreport m_EFieldModRepx, m_EFieldModRepy, m_EFieldModRepz, m_EFieldPhaseRepx, m_EFieldPhaseRepy, m_EFieldPhaseRepz, m_BFieldModRepx, m_BFieldModRepy, m_BFieldModRepz, m_BFieldPhaseRepx, m_BFieldPhaseRepy, m_BFieldPhaseRepz;
    alglib::real_2d_array m_Points_and_EModx, m_Points_and_EMody, m_Points_and_EModz, m_Points_and_EPhasex, m_Points_and_EPhasey, m_Points_and_EPhasez, m_Points_and_BModx, m_Points_and_BMody, m_Points_and_BModz, m_Points_and_BPhasex, m_Points_and_BPhasey, m_Points_and_BPhasez;
    // std::array<std::array<double, 6>, n_points> ;
};

template <const size_t n_points>
class FieldInterpolatedRIComps : public FieldBase
{
public:
    FieldInterpolatedRIComps(double lambda, std::array<double, n_points> &xs, std::array<double, n_points> &ys, std::array<double, n_points> &zs, double PointsDist, std::array<FieldValue, n_points> &FieldInRegion) :
        FieldBase(lambda), m_npoints(n_points), m_xs(xs),  m_ys(ys),  m_zs(zs), m_FieldInRegion(FieldInRegion), m_PointsDist(PointsDist)

    {
        alglib::rbfcreate(3,1, m_EFieldRealx);
        alglib::rbfcreate(3,1, m_EFieldRealy);
        alglib::rbfcreate(3,1, m_EFieldRealz);
        alglib::rbfcreate(3,1, m_BFieldRealx);
        alglib::rbfcreate(3,1, m_BFieldRealy);
        alglib::rbfcreate(3,1, m_BFieldRealz);
        alglib::rbfcreate(3,1, m_EFieldImagx);
        alglib::rbfcreate(3,1, m_EFieldImagy);
        alglib::rbfcreate(3,1, m_EFieldImagz);
        alglib::rbfcreate(3,1, m_BFieldImagx);
        alglib::rbfcreate(3,1, m_BFieldImagy);
        alglib::rbfcreate(3,1, m_BFieldImagz);

        m_Points_and_ERealx.setlength(m_npoints, 4);
        m_Points_and_ERealy.setlength(m_npoints, 4);
        m_Points_and_ERealz.setlength(m_npoints, 4);
        m_Points_and_EImagx.setlength(m_npoints, 4);
        m_Points_and_EImagy.setlength(m_npoints, 4);
        m_Points_and_EImagz.setlength(m_npoints, 4);
        m_Points_and_BRealx.setlength(m_npoints, 4);
        m_Points_and_BRealy.setlength(m_npoints, 4);
        m_Points_and_BRealz.setlength(m_npoints, 4);
        m_Points_and_BImagx.setlength(m_npoints, 4);
        m_Points_and_BImagy.setlength(m_npoints, 4);
        m_Points_and_BImagz.setlength(m_npoints, 4);

       
        for (int i=0; i<m_npoints; i++) {
            m_Points_and_ERealx[i][0]=m_xs[i];
            m_Points_and_ERealx[i][1]=m_ys[i];
            m_Points_and_ERealx[i][2]=m_zs[i];
            m_Points_and_ERealx[i][3]=vec_real(m_FieldInRegion[i].E)[0];

            m_Points_and_ERealy[i][0]=m_xs[i];
            m_Points_and_ERealy[i][1]=m_ys[i];
            m_Points_and_ERealy[i][2]=m_zs[i];
            m_Points_and_ERealy[i][3]=vec_real(m_FieldInRegion[i].E)[1];

            m_Points_and_ERealz[i][0]=m_xs[i];
            m_Points_and_ERealz[i][1]=m_ys[i];
            m_Points_and_ERealz[i][2]=m_zs[i];
            m_Points_and_ERealz[i][3]=vec_real(m_FieldInRegion[i].E)[2];

            m_Points_and_EImagx[i][0]=m_xs[i];
            m_Points_and_EImagx[i][1]=m_ys[i];
            m_Points_and_EImagx[i][2]=m_zs[i];
            m_Points_and_EImagx[i][3]=vec_imag(m_FieldInRegion[i].E)[0];

            m_Points_and_EImagy[i][0]=m_xs[i];
            m_Points_and_EImagy[i][1]=m_ys[i];
            m_Points_and_EImagy[i][2]=m_zs[i];
            m_Points_and_EImagy[i][3]=vec_imag(m_FieldInRegion[i].E)[1];

            m_Points_and_EImagz[i][0]=m_xs[i];
            m_Points_and_EImagz[i][1]=m_ys[i];
            m_Points_and_EImagz[i][2]=m_zs[i];
            m_Points_and_EImagz[i][3]=vec_imag(m_FieldInRegion[i].E)[2];

            m_Points_and_BRealx[i][0]=m_xs[i];
            m_Points_and_BRealx[i][1]=m_ys[i];
            m_Points_and_BRealx[i][2]=m_zs[i];
            m_Points_and_BRealx[i][3]=vec_real(m_FieldInRegion[i].B)[0];

            m_Points_and_BRealy[i][0]=m_xs[i];
            m_Points_and_BRealy[i][1]=m_ys[i];
            m_Points_and_BRealy[i][2]=m_zs[i];
            m_Points_and_BRealy[i][3]=vec_real(m_FieldInRegion[i].B)[1];

            m_Points_and_BRealz[i][0]=m_xs[i];
            m_Points_and_BRealz[i][1]=m_ys[i];
            m_Points_and_BRealz[i][2]=m_zs[i];
            m_Points_and_BRealz[i][3]=vec_real(m_FieldInRegion[i].B)[2];

            m_Points_and_BImagx[i][0]=m_xs[i];
            m_Points_and_BImagx[i][1]=m_ys[i];
            m_Points_and_BImagx[i][2]=m_zs[i];
            m_Points_and_BImagx[i][3]=vec_imag(m_FieldInRegion[i].E)[0];

            m_Points_and_BImagy[i][0]=m_xs[i];
            m_Points_and_BImagy[i][1]=m_ys[i];
            m_Points_and_BImagy[i][2]=m_zs[i];
            m_Points_and_BImagy[i][3]=vec_imag(m_FieldInRegion[i].E)[1];

            m_Points_and_BImagz[i][0]=m_xs[i];
            m_Points_and_BImagz[i][1]=m_ys[i];
            m_Points_and_BImagz[i][2]=m_zs[i];
            m_Points_and_BImagz[i][3]=vec_imag(m_FieldInRegion[i].E)[2];
        }

        alglib::rbfsetpoints(m_EFieldRealx, m_Points_and_ERealx);
        alglib::rbfsetpoints(m_EFieldRealy, m_Points_and_ERealy);
        alglib::rbfsetpoints(m_EFieldRealz, m_Points_and_ERealz);
        alglib::rbfsetpoints(m_EFieldImagx, m_Points_and_EImagx);
        alglib::rbfsetpoints(m_EFieldImagy, m_Points_and_EImagy);
        alglib::rbfsetpoints(m_EFieldImagz, m_Points_and_EImagz);
        alglib::rbfsetpoints(m_BFieldRealx, m_Points_and_BRealx);
        alglib::rbfsetpoints(m_BFieldRealy, m_Points_and_BRealy);
        alglib::rbfsetpoints(m_BFieldRealz, m_Points_and_BRealz);
        alglib::rbfsetpoints(m_BFieldImagx, m_Points_and_BImagx);
        alglib::rbfsetpoints(m_BFieldImagy, m_Points_and_BImagy);
        alglib::rbfsetpoints(m_BFieldImagz, m_Points_and_BImagz);


        alglib::rbfsetalgohierarchical(m_EFieldRealx, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_EFieldRealx, 0.0);
        alglib::rbfbuildmodel(m_EFieldRealx, m_EFieldRealRepx);

        alglib::rbfsetalgohierarchical(m_EFieldRealy, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_EFieldRealy, 0.0);
        alglib::rbfbuildmodel(m_EFieldRealy, m_EFieldRealRepy);

        alglib::rbfsetalgohierarchical(m_EFieldRealz, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_EFieldRealz, 0.0);
        alglib::rbfbuildmodel(m_EFieldRealz, m_EFieldRealRepz);

        alglib::rbfsetalgohierarchical(m_EFieldImagx, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_EFieldImagx,0.0);
        alglib::rbfbuildmodel(m_EFieldImagx, m_EFieldImagRepx);

        alglib::rbfsetalgohierarchical(m_EFieldImagy, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_EFieldImagy,0.0);
        alglib::rbfbuildmodel(m_EFieldImagy, m_EFieldImagRepy);

        alglib::rbfsetalgohierarchical(m_EFieldImagz, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_EFieldImagz,0.0);
        alglib::rbfbuildmodel(m_EFieldImagz, m_EFieldImagRepz);

        alglib::rbfsetalgohierarchical(m_BFieldRealx, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_BFieldRealx, 0.0);
        alglib::rbfbuildmodel(m_BFieldRealx, m_BFieldRealRepx);

        alglib::rbfsetalgohierarchical(m_BFieldRealy, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_BFieldRealy, 0.0);
        alglib::rbfbuildmodel(m_BFieldRealy, m_BFieldRealRepy);

        alglib::rbfsetalgohierarchical(m_BFieldRealz, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_BFieldRealz, 0.0);
        alglib::rbfbuildmodel(m_BFieldRealz, m_BFieldRealRepz);

        alglib::rbfsetalgohierarchical(m_BFieldImagx, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_BFieldImagx, 0.0);
        alglib::rbfbuildmodel(m_BFieldImagx, m_BFieldImagRepx);

        alglib::rbfsetalgohierarchical(m_BFieldImagy, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_BFieldImagy, 0.0);
        alglib::rbfbuildmodel(m_BFieldImagy, m_BFieldImagRepy);

        alglib::rbfsetalgohierarchical(m_BFieldImagz, m_PointsDist, 7, 0.0);
        // alglib::rbfsetalgothinplatespline(m_BFieldImagz, 0.0);
        alglib::rbfbuildmodel(m_BFieldImagz, m_BFieldImagRepz);
        

        std::cout << "ERealx check (should be 1) ---> " << int(m_EFieldRealRepx.terminationtype) << std::endl;
        std::cout << "ERealy check (should be 1) ---> " << int(m_EFieldRealRepy.terminationtype) << std::endl;
        std::cout << "ERealz check (should be 1) ---> " << int(m_EFieldRealRepz.terminationtype) << std::endl;
        std::cout << "EImagx check (should be 1) ---> " << int(m_EFieldImagRepx.terminationtype) << std::endl;
        std::cout << "EImagy check (should be 1) ---> " << int(m_EFieldImagRepy.terminationtype) << std::endl;
        std::cout << "EImagz check (should be 1) ---> " << int(m_EFieldImagRepz.terminationtype) << std::endl;
        std::cout << "BRealx check (should be 1) ---> " << int(m_BFieldRealRepx.terminationtype) << std::endl;
        std::cout << "BRealy check (should be 1) ---> " << int(m_BFieldRealRepy.terminationtype) << std::endl;
        std::cout << "BRealz check (should be 1) ---> " << int(m_BFieldRealRepz.terminationtype) << std::endl;
        std::cout << "BImagx check (should be 1) ---> " << int(m_BFieldImagRepx.terminationtype) << std::endl;
        std::cout << "BImagy check (should be 1) ---> " << int(m_BFieldImagRepy.terminationtype) << std::endl;
        std::cout << "BImagz check (should be 1) ---> " << int(m_BFieldImagRepz.terminationtype) << std::endl;
    }

    
    FieldValue get(const Position& pos) const override
    {
        FieldValue result;
        // alglib::real_1d_array coord, EReal, EImag, BReal, BImag;
        alglib::real_1d_array coord, ERealx, ERealy, ERealz, EImagx, EImagy, EImagz, BRealx, BRealy, BRealz, BImagx, BImagy, BImagz;
        coord.setlength(3);
        ERealx.setlength(1);
        ERealy.setlength(1);
        ERealz.setlength(1);
        EImagx.setlength(1);
        EImagy.setlength(1);
        EImagz.setlength(1);
        BRealx.setlength(1);
        BRealy.setlength(1);
        BRealz.setlength(1);
        BImagx.setlength(1);
        BImagy.setlength(1);
        BImagz.setlength(1);

        coord[0]=pos[0];
        coord[1]=pos[1];
        coord[2]=pos[2];

        alglib::rbfcalcbuffer ERBUFx, ERBUFy, ERBUFz, EIBUFx, EIBUFy, EIBUFz, BRBUFx, BRBUFy, BRBUFz, BIBUFx, BIBUFy, BIBUFz;

        alglib::rbfcreatecalcbuffer(m_EFieldRealx, ERBUFx);
        alglib::rbfcreatecalcbuffer(m_EFieldRealy, ERBUFy);
        alglib::rbfcreatecalcbuffer(m_EFieldRealz, ERBUFz);
        alglib::rbfcreatecalcbuffer(m_EFieldImagx, EIBUFx);
        alglib::rbfcreatecalcbuffer(m_EFieldImagy, EIBUFy);
        alglib::rbfcreatecalcbuffer(m_EFieldImagz, EIBUFz);
        alglib::rbfcreatecalcbuffer(m_BFieldRealx, BRBUFx);
        alglib::rbfcreatecalcbuffer(m_BFieldRealy, BRBUFy);
        alglib::rbfcreatecalcbuffer(m_BFieldRealz, BRBUFz);
        alglib::rbfcreatecalcbuffer(m_BFieldImagx, BIBUFx);
        alglib::rbfcreatecalcbuffer(m_BFieldImagy, BIBUFy);
        alglib::rbfcreatecalcbuffer(m_BFieldImagz, BIBUFz);
        
        alglib::rbftscalcbuf(m_EFieldRealx, ERBUFx, coord, ERealx);
        alglib::rbftscalcbuf(m_EFieldRealy, ERBUFy, coord, ERealy);
        alglib::rbftscalcbuf(m_EFieldRealz, ERBUFz, coord, ERealz);
        alglib::rbftscalcbuf(m_EFieldImagx, EIBUFx, coord, EImagx);
        alglib::rbftscalcbuf(m_EFieldImagy, EIBUFy, coord, EImagy);
        alglib::rbftscalcbuf(m_EFieldImagz, EIBUFz, coord, EImagz);
        alglib::rbftscalcbuf(m_BFieldRealx, BRBUFx, coord, BRealx);
        alglib::rbftscalcbuf(m_BFieldRealy, BRBUFy, coord, BRealy);
        alglib::rbftscalcbuf(m_BFieldRealz, BRBUFz, coord, BRealz);
        alglib::rbftscalcbuf(m_BFieldImagx, BIBUFx, coord, BImagx);
        alglib::rbftscalcbuf(m_BFieldImagy, BIBUFy, coord, BImagy);
        alglib::rbftscalcbuf(m_BFieldImagz, BIBUFz, coord, BImagz);

        result.E[0] = ERealx[0]+std::complex<double>(0.0, 1.0)*EImagx[0];
        result.E[1] = ERealy[0]+std::complex<double>(0.0, 1.0)*EImagy[0];
        result.E[2] = ERealz[0]+(std::complex<double>(0.0, 1.0)*EImagz[0]);
        result.B[0] = BRealx[0]+(std::complex<double>(0.0, 1.0)*BImagx[0]);
        result.B[1] = BRealy[0]+(std::complex<double>(0.0, 1.0)*BImagy[0]);
        result.B[2] = BRealz[0]+(std::complex<double>(0.0, 1.0)*BImagz[0]);

        return result;
    }

private:
    size_t m_npoints;
    double m_PointsDist;
    std::array<double, n_points> m_xs, m_ys, m_zs;
    std::array<FieldValue, n_points> &m_FieldInRegion;
    alglib::rbfmodel m_EFieldRealx, m_EFieldRealy, m_EFieldRealz, m_EFieldImagx, m_EFieldImagy, m_EFieldImagz, m_BFieldRealx, m_BFieldRealy, m_BFieldRealz, m_BFieldImagx, m_BFieldImagy, m_BFieldImagz;
    alglib::rbfreport m_EFieldRealRepx, m_EFieldRealRepy, m_EFieldRealRepz, m_EFieldImagRepx, m_EFieldImagRepy, m_EFieldImagRepz, m_BFieldRealRepx, m_BFieldRealRepy, m_BFieldRealRepz, m_BFieldImagRepx, m_BFieldImagRepy, m_BFieldImagRepz;
    alglib::real_2d_array m_Points_and_ERealx, m_Points_and_ERealy, m_Points_and_ERealz, m_Points_and_EImagx, m_Points_and_EImagy, m_Points_and_EImagz, m_Points_and_BRealx, m_Points_and_BRealy, m_Points_and_BRealz, m_Points_and_BImagx, m_Points_and_BImagy, m_Points_and_BImagz;
    // std::array<std::array<double, 6>, n_points> ;
};


*/

#endif // FIELD_HPP
