#include "stratton-chu/csv-saver.hpp"

CSVSaver::CSVSaver(const char* filename) :
    m_ofs(filename)
{
    m_ofs << "x,y,z,Ex_re,Ex_im,Ey_re,Ey_im,Ez_re,Ez_im,Bx_re,Bx_im,By_re,By_im,Bz_re,Bz_im" << std::endl;
}

void CSVSaver::add_point(Position p, FieldValue fv)
{
    m_ofs << p[0] << "," << p[1] << ","<< p[2] << ","

          << fv.E[0].real() << "," << fv.E[0].imag() << ","
          << fv.E[1].real() << "," << fv.E[1].imag() << ","
          << fv.E[2].real() << "," << fv.E[2].imag() << ","

          << fv.B[0].real() << "," << fv.B[0].imag() << ","
          << fv.B[1].real() << "," << fv.B[1].imag() << ","
          << fv.B[2].real() << "," << fv.B[2].imag() << std::endl << std::flush;
}
