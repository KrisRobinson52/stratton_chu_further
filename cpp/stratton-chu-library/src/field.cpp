#include "stratton-chu/field.hpp"
#include "stratton-chu/geometry.hpp"

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

// template <const size_t n_points>
// FieldInterpolated::FieldInterpolated (double lambda, std::vector<double> &xs, std::vector<double> &ys, std::vector<double> &zs, double PointsDist, std::vector<FieldValue> &FieldInRegion) :
//     FieldBase(lambda), m_xs(xs),  m_ys(ys),  m_zs(zs), m_FieldInRegion(FieldInRegion), m_PointsDist(PointsDist)
// {
//     alglib::rbfcreate(3,3, m_EFieldReal);
//     alglib::rbfcreate(3,3, m_EFieldImag);
//     alglib::rbfcreate(3,3, m_BFieldReal);
//     alglib::rbfcreate(3,3, m_BFieldImag);
//     points_numb=xs.size();
//     for (int i=0; i<points_numb; i++) {
//         m_Points_and_EReals[i][0]=m_xs[i];
//         m_Points_and_EReals[i][1]=m_ys[i];
//         m_Points_and_EReals[i][2]=m_zs[i];
//         m_Points_and_EReals[i][3]=vec_real(m_FieldInRegion[i].E)[0];
//         m_Points_and_EReals[i][4]=vec_real(m_FieldInRegion[i].E)[1];
//         m_Points_and_EReals[i][5]=vec_real(m_FieldInRegion[i].E)[2];

//         m_Points_and_EImags[i][0]=m_xs[i];
//         m_Points_and_EImags[i][1]=m_ys[i];
//         m_Points_and_EImags[i][2]=m_zs[i];
//         m_Points_and_EImags[i][3]=vec_imag(m_FieldInRegion[i].E)[0];
//         m_Points_and_EImags[i][4]=vec_imag(m_FieldInRegion[i].E)[1];
//         m_Points_and_EImags[i][5]=vec_imag(m_FieldInRegion[i].E)[2];

//         m_Points_and_BReals[i][0]=m_xs[i];
//         m_Points_and_BReals[i][1]=m_ys[i];
//         m_Points_and_BReals[i][2]=m_zs[i];
//         m_Points_and_BReals[i][3]=vec_real(m_FieldInRegion[i].B)[0];
//         m_Points_and_BReals[i][4]=vec_real(m_FieldInRegion[i].B)[1];
//         m_Points_and_BReals[i][5]=vec_real(m_FieldInRegion[i].B)[2];

//         m_Points_and_BImags[i][0]=m_xs[i];
//         m_Points_and_BImags[i][1]=m_ys[i];
//         m_Points_and_BImags[i][2]=m_zs[i];
//         m_Points_and_BImags[i][3]=vec_imag(m_FieldInRegion[i].B)[0];
//         m_Points_and_BImags[i][4]=vec_imag(m_FieldInRegion[i].B)[1];
//         m_Points_and_BImags[i][5]=vec_imag(m_FieldInRegion[i].B)[2];
//     }
//     alglib::rbfsetpoints(m_EFieldReal, m_Points_and_EReals);
//     alglib::rbfsetpoints(m_EFieldImag, m_Points_and_EImags);
//     alglib::rbfsetpoints(m_BFieldReal, m_Points_and_BReals);
//     alglib::rbfsetpoints(m_BFieldImag, m_Points_and_BImags);

//     alglib::rbfsetalgohierarchical(m_EFieldReal, m_PointsDist, 1, 0.0);
//     alglib::rbfbuildmodel(m_EFieldReal, m_EFieldRealRep);

//     alglib::rbfsetalgohierarchical(m_EFieldImag, m_PointsDist, 1, 0.0);
//     alglib::rbfbuildmodel(m_EFieldImag, m_EFieldImagRep);

//     alglib::rbfsetalgohierarchical(m_BFieldReal, m_PointsDist, 1, 0.0);
//     alglib::rbfbuildmodel(m_BFieldReal, m_BFieldRealRep);

//     alglib::rbfsetalgohierarchical(m_BFieldImag, m_PointsDist, 1, 0.0);
//     alglib::rbfbuildmodel(m_BFieldImag, m_BFieldImagRep);
    

//     std::cout << "EReal check (should be 1) ---> " << int(m_EFieldRealRep.terminationtype) << std::endl;
//     std::cout << "EImag check (should be 1) ---> " << int(m_EFieldImagRep.terminationtype) << std::endl;
//     std::cout << "BReal check (should be 1) ---> " << int(m_BFieldRealRep.terminationtype) << std::endl;
//     std::cout << "BImag check (should be 1) ---> " << int(m_BFieldImagRep.terminationtype) << std::endl;    
// }

// FieldValue FieldInterpolated::get(const Position& pos) const
// {
//     FieldValue result;
//     std::array<double, 3> coord;
//     coord[0]=pos[0];
//     coord[1]=pos[1];
//     coord[2]=pos[2];

//     std::array<double, 3> EReal, EImag, BReal, BImag;
//     alglib::rbfcalc(m_EFieldReal, coord, EReal);
//     alglib::rbfcalc(m_EFieldImag, coord, EImag);
//     alglib::rbfcalc(m_BFieldReal, coord, BReal);
//     alglib::rbfcalc(m_BFieldImag, coord, BImag);


//     result.E[0] = EReal[0]+std::complex<double>(0.0, 1.0)*EImag[0];
//     result.E[1] = EReal[1]+std::complex<double>(0.0, 1.0)*EImag[1];
//     result.E[2] = EReal[1]+std::complex<double>(0.0, 1.0)*EImag[2];
//     result.B[0] = BReal[0]+std::complex<double>(0.0, 1.0)*BImag[0];
//     result.B[1] = BReal[1]+std::complex<double>(0.0, 1.0)*BImag[1];
//     result.B[2] = BReal[2]+std::complex<double>(0.0, 1.0)*BImag[2];

//     return result;
// }