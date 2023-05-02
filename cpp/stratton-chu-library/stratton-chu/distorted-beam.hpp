#ifndef DISTORTEDBEAM_HPP
#define DISTORTEDBEAM_HPP

#include <vector>
#include <array>
#include <random>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

#include <zernike_bits/zernike_radial_poly.h>
#include <zernike_bits/zernike_poly.h>


double ZernikeCartesian(unsigned int n, int m, double x1, double x2, double beam_rad){
    double r = sqrt(x1*x1+x2*x2)/beam_rad;
    double theta=0;
    if (r!=0.0){
        theta=atan2(x2 ,x1);
    }   
    return Zernike::ZernikePolynomial(n,m,r,theta);
}

template<size_t Order>
struct ZernikeAberrations
{
    static std::array<double, (Order + 1) * (Order + 1)> gen_ampls(double sigma, size_t beam_ind);

    const std::array<double, (Order + 1) * (Order + 1)> m_ampls;
    const double m_beam_radius;

    ZernikeAberrations(double sigma, double beam_radius, size_t beam_ind)
        : m_ampls{gen_ampls(sigma, beam_ind)}, m_beam_radius{beam_radius}
    {}

    double operator()(double x1, double x2) const;
};

template<size_t Order>
double ZernikeAberrations<Order>::operator()(double x1, double x2) const
{
    std::vector<double> m_aberr;
    for (long n = Order; n > 0; n--) {
        for (long m = -n; m <= n; m++) {
            if ((n - m) % 2 == 0) {
                m_aberr.push_back(ZernikeCartesian(n, m, x1, x2, m_beam_radius));
            }
        }
    }

    size_t aberr_ind = 0;
    for (size_t i = Order; i > 0; i--) {
        for (size_t j = 0; j < i + 1; j++) {
            double ampl_in_mode = m_ampls[i * Order + j];
            m_aberr[aberr_ind++] *= ampl_in_mode;
        }
    }

    double summ=0;
    for (size_t i=0; i<m_aberr.size(); i++){
        summ+=m_aberr[i];
    }
    return summ;
    // return std::accumulate(m_aberr.begin(), m_aberr.end(), 0);
}

template<size_t Order>
std::array<double, (Order + 1) * (Order + 1)> ZernikeAberrations<Order>::gen_ampls(double sigma, size_t beam_ind)
{
    std::string fname="poly_ampl"+std::to_string(beam_ind)+".txt";
    std::ofstream outr(fname);
    std::random_device dev;
    std::mt19937 rng(dev());
    std::array<double, (Order + 1) * (Order + 1)> m_ampls;
    std::normal_distribution<double> dist{0, sigma};
    size_t mode_num=0;
    double ampls_sum=0.0;
    std::vector<double> m_diffs;
    for (size_t i = Order; i > 0; i--) {
        for (size_t j = 0; j < i + 1; j++) {
            double ampl_in_mode = dist(rng);
            ampls_sum+=ampl_in_mode;
            m_ampls[i * Order + j] = ampl_in_mode;
            outr << "mode " << i << " " << j << " with coef "
                 << ampl_in_mode << '\n';
            mode_num++;
        }
    }
    double avg=ampls_sum/mode_num;
    double diffs_sum=0.0;
    for (size_t i = Order; i > 0; i--) {
        for (size_t j = 0; j < i + 1; j++) {
            m_diffs.push_back(sqr(m_ampls[i * Order + j]-avg));
            diffs_sum+=sqr(m_ampls[i * Order + j]-avg);
        }
    }
    double sko=sqrt(diffs_sum/mode_num);
    outr << sko/M_PI*180 << '\n';
    // for (size_t n = 0; n < (Order + 1) * (Order + 1); n++) {
    //     outr << m_ampls[n] << '\n';
    // }

    // for (long n = Order; n > 0; n--) {
    //     for (long m = -n; m <= n; m++) {
    //         if ((n - m) % 2 == 0) {
    //             outr << "mode " << n << " " << m << " says hi"<< '\n';
    //         }
    //     }
    // }
    return m_ampls;
}



#endif // DISTORTEDBEAM_HPP