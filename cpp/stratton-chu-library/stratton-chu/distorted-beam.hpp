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
    double norm = std::sqrt((2.0*n+2.0)/((1.0+(m==0 ? 1.0 : 0.0))));
    // if (r!=0.0){
    if (r>1e-20){
        theta=atan2(x2 ,x1);
    }
    return Zernike::ZernikePolynomial(n,m,r,theta);   
    // return Zernike::ZernikePolynomial(n,m,r,theta)/norm;
}

double ZernikeCartesianOnSquare1(unsigned int j, double x1, double x2, double beam_rad){
    double newx1 = x1/beam_rad/sqrt(2);
    double newx2 = x2/beam_rad/sqrt(2);
    
    double rsqr = newx1*newx1+newx2*newx2;

    double ZerPol;

    switch(j){
        case 1:
            ZerPol=1;
            break;
        case 2:
            ZerPol=sqrt(6)*newx1;
            break;
        case 3:
            ZerPol=sqrt(6)*newx2;
            break;
        case 4:
            ZerPol=sqrt(5.0/2.0)*(3*rsqr-1);
            break;
        case 5:
            ZerPol=6*newx1*newx2;
            break;
        case 6:
            ZerPol=3*sqrt(5.0/2.0)*(newx1*newx1-newx2*newx2);
            break;
        case 7:
            ZerPol=sqrt(21.0/31.0)*(15*rsqr-7)*newx2;
            break;
        case 8:
            ZerPol=sqrt(21.0/31.0)*(15*rsqr-7)*newx1;
            break;
        case 9:
            ZerPol=sqrt(5.0/31.0)*(27*newx1*newx1-35*newx2*newx2+6)*newx2;
            break;
        case 10:
            ZerPol=sqrt(5.0/31.0)*(35*newx1*newx1-27*newx2*newx2-6)*newx1;
            break;
        case 11:
            ZerPol=1.0/(2*sqrt(67))*(315*rsqr*rsqr-240*rsqr+31);
            break;
        case 12:
            ZerPol=15.0/(2*sqrt(2))*(newx1*newx1-newx2*newx2)*(7*rsqr-3);
            break;
        case 13:
            ZerPol=sqrt(42)*(5*rsqr-3)*newx1*newx2;
            break;
        case 14:
            ZerPol=3.0/(4*sqrt(134))*(10*(49*pow(newx1, 4)-36*newx1*newx1*newx2*newx2+49*pow(newx2, 4))-150*rsqr+11);
            break;
        case 15:
            ZerPol=5*sqrt(42)*(newx1*newx1-newx2*newx2)*newx1*newx2;
            break;
        case 16:
            ZerPol=sqrt(55.0/1966.0)*(315*rsqr*rsqr-280*newx1*newx1-324*newx2*newx2+57)*newx1;
            break;
        case 17:
            ZerPol=sqrt(55.0/1966.0)*(315*rsqr*rsqr-324*newx1*newx1-280*newx2*newx2+57)*newx2;
            break;
        case 18:
            ZerPol=0.5*sqrt(3.0/844397.0)*(105*(1023*pow(newx1, 4)+80*newx1*newx1*newx2*newx2-943*pow(newx2, 4))-61075*newx1*newx1+39915*newx2*newx2+4692)*newx1;
            break;
        case 19:
            ZerPol=0.5*sqrt(3.0/844397.0)*(105*(943*pow(newx1, 4)-80*newx1*newx1*newx2*newx2-1023*pow(newx2, 4))-39915*newx1*newx1+61075*newx2*newx2-4692)*newx2;
            break;
        case 20:
            ZerPol=0.25*sqrt(7.0/859.0)*(6*(693*pow(newx1, 4)-500*newx1*newx1*newx2*newx2+525*pow(newx2, 4))-1810*newx1*newx1-450*newx2*newx2+165)*newx1;
            break;
        case 21:
            ZerPol=0.25*sqrt(7.0/859.0)*(6*(525*pow(newx1, 4)-500*newx1*newx1*newx2*newx2+693*pow(newx2, 4))-450*newx1*newx1-1810*newx2*newx2+165)*newx2;
            break;
        case 22:
            ZerPol=0.25*sqrt(65.0/849.0)*(1155*rsqr*rsqr*rsqr-15*(91*pow(newx1, 4)+198*newx1*newx1*newx2*newx2+91*pow(newx2, 4))+453*rsqr-31);
            break;
        case 23:
            ZerPol=sqrt(33.0/3923.0)*(1575*rsqr*rsqr-1820*rsqr+471)*newx1*newx2;
            break;
        case 24:
            ZerPol=21.0/4.0*sqrt(65.0/1349.0)*(165*rsqr*rsqr-140*rsqr+27)*(newx1*newx1-newx2*newx2);
            break;
        case 25:
            ZerPol=7*sqrt(33.0/2.0)*(9*rsqr-5)*newx1*newx2*(newx1*newx1-newx2*newx2);
            break;
        case 26:
            ZerPol=1.0/8.0/sqrt(849.0)*(42*(1573*pow(newx1, 6)-375*pow(newx1, 4)*newx2*newx2-375*pow(newx2, 4)*newx1*newx1+1573*pow(newx2, 6))-60*(707*pow(newx1, 4)-225*newx1*newx1*newx2*newx2+707*pow(newx2, 4))+6045*rsqr-245);
            break;
        case 27:
            ZerPol=0.5/sqrt(7846)*(14*(2673*pow(newx1, 4)-2500*newx1*newx1*newx2*newx2+2673*pow(newx2, 4))-10290*rsqr+1305)*newx1*newx2;
            break;
        case 28:
            ZerPol=21.0/8.0/sqrt(1349.0)*(3146*pow(newx1, 6)-2250*pow(newx1, 4)*newx2*newx2-2250*pow(newx2, 4)*newx1*newx1-3146*pow(newx2, 6)-1770*(pow(newx1, 4)-pow(newx2, 4))+245*(newx1*newx1-newx2*newx2));
            break;
        default:
            throw std::runtime_error("Unidentified Zernike index");
    }

    
    return ZerPol;   
    // return Zernike::ZernikePolynomial(n,m,r,theta)/norm;
}

double ZernikeCartesianOnSquare2(unsigned int j, double x1, double x2, double beam_rad){
    double newx1 = x1/beam_rad;
    double newx2 = x2/beam_rad;
    
    double rsqr = newx1*newx1+newx2*newx2;

    double ZerPol;

    switch(j){
        case 1:
            ZerPol=1;
            break;
        case 2:
            ZerPol=sqrt(3)*newx1;
            break;
        case 3:
            ZerPol=sqrt(3)*newx2;
            break;
        case 4:
            ZerPol=0.5*sqrt(5.0/2.0)*(3*rsqr-2);
            break;
        case 5:
            ZerPol=3*newx1*newx2;
            break;
        case 6:
            ZerPol=1.5*sqrt(5.0/2.0)*(newx1*newx1-newx2*newx2);
            break;
        case 7:
            ZerPol=0.5*sqrt(21.0/62.0)*(15*rsqr-14)*newx2;
            break;
        case 8:
            ZerPol=0.5*sqrt(21.0/62.0)*(15*rsqr-14)*newx1;
            break;
        case 9:
            ZerPol=0.5*sqrt(5.0/62.0)*(27*newx1*newx1-35*newx2*newx2+12)*newx2;
            break;
        case 10:
            ZerPol=0.5*sqrt(5.0/62.0)*(35*newx1*newx1-27*newx2*newx2-12)*newx1;
            break;
        case 11:
            ZerPol=1.0/(8*sqrt(67))*(315*rsqr*rsqr-480*rsqr+124);
            break;
        case 12:
            ZerPol=15.0/(8*sqrt(2))*(newx1*newx1-newx2*newx2)*(7*rsqr-6);
            break;
        case 13:
            ZerPol=0.5*sqrt(21.0/2.0)*(5*rsqr-6)*newx1*newx2;
            break;
        case 14:
            ZerPol=3.0/(8*sqrt(134))*(5*(49*pow(newx1, 4)-36*newx1*newx1*newx2*newx2+49*pow(newx2, 4))-150*rsqr+22);
            break;
        case 15:
            ZerPol=2.5*sqrt(21.0/2.0)*(newx1*newx1-newx2*newx2)*newx1*newx2;
            break;
        case 16:
            ZerPol=0.25*sqrt(55.0/3932.0)*(315*rsqr*rsqr-560*newx1*newx1-648*newx2*newx2+228)*newx1;
            break;
        case 17:
            ZerPol=0.25*sqrt(55.0/3932.0)*(315*rsqr*rsqr-648*newx1*newx1-560*newx2*newx2+228)*newx2;
            break;
        case 18:
            ZerPol=0.125*sqrt(3.0/1688794.0)*(105*(1023*pow(newx1, 4)+80*newx1*newx1*newx2*newx2-943*pow(newx2, 4))-122150*newx1*newx1+79830*newx2*newx2+18768)*newx1;
            break;
        case 19:
            ZerPol=0.125*sqrt(3.0/1688794.0)*(105*(943*pow(newx1, 4)-80*newx1*newx1*newx2*newx2-1023*pow(newx2, 4))-79830*newx1*newx1+122150*newx2*newx2-18768)*newx2;
            break;
        case 20:
            ZerPol=0.125*sqrt(7.0/1718.0)*(3*(693*pow(newx1, 4)-500*newx1*newx1*newx2*newx2+525*pow(newx2, 4))-1810*newx1*newx1-450*newx2*newx2+330)*newx1;
            break;
        case 21:
            ZerPol=0.125*sqrt(7.0/1718.0)*(3*(525*pow(newx1, 4)-500*newx1*newx1*newx2*newx2+693*pow(newx2, 4))-450*newx1*newx1-1810*newx2*newx2+330)*newx2;
            break;
        // case 16:
        //     ZerPol=(9.313762357*rsqr*rsqr-16.55779969*newx1*newx1-19.15973967*newx2*newx2+6.741389845)*newx1;
        //     break;
        // case 17:
        //     ZerPol=(9.313762357*rsqr*rsqr-19.15973967*newx1*newx1-16.55779969*newx2*newx2+6.741389845)*newx2;
        //     break;
        // case 18:
        //     ZerPol=(17.89564012*pow(newx1, 4)+1.399463791*newx1*newx1*newx2*newx2-16.49617633*pow(newx2, 4)-20.35053247*newx1*newx1+13.29990146*newx2*newx2+3.126801448)*newx1;
        //     break;
        // case 19:
        //     ZerPol=(16.49617633*pow(newx1, 4)-1.399463791*newx1*newx1*newx2*newx2-17.89564012*pow(newx2, 4)-13.29990146*newx1*newx1+20.35053247*newx2*newx2-3.126801448)*newx2;
        //     break;
        // case 20:
        //     ZerPol=(16.58830278*pow(newx1, 4)-11.96847179*newx1*newx1*newx2*newx2+12.56689562*pow(newx2, 4)-14.44195688*newx1*newx1-3.590541684*newx2*newx2+2.633063969)*newx1;
        //     break;
        // case 21:
        //     ZerPol=(12.56689562*pow(newx1, 4)-11.96847179*newx1*newx1*newx2*newx2+16.58830278*pow(newx2, 4)-3.590541684*newx1*newx1-14.44195688*newx2*newx2+2.633063969)*newx2;
        //     break;
        default:
            throw std::runtime_error("Unidentified Zernike index");
    }

    
    return ZerPol;   
    // return Zernike::ZernikePolynomial(n,m,r,theta)/norm;
}

double ZernikeCartesianOnSquare3(unsigned int j, double x1, double x2, double beam_rad){
    double newx1 = x1/beam_rad*sqrt(M_PI)/2.0;
    double newx2 = x2/beam_rad*sqrt(M_PI)/2.0;
    
    double rsqr = newx1*newx1+newx2*newx2;

    double ZerPol;

    switch(j){
        case 1:
            ZerPol=1;
            break;
        case 2:
            ZerPol=1.95*newx1;
            break;
        case 3:
            ZerPol=1.95*newx2;
            break;
        case 4:
            ZerPol=3.02*rsqr-1.59;
            break;
        case 5:
            ZerPol=3.82*newx1*newx2;
            break;
        case 6:
            ZerPol=3.03*(newx1*newx1-newx2*newx2);
            break;
        case 7:
            ZerPol=(6.24*rsqr-4.58)*newx2;
            break;
        case 8:
            ZerPol=(6.24*rsqr-4.58)*newx1;
            break;
        case 9:
            ZerPol=(5.5*newx1*newx1-7.06*newx2*newx2+1.9)*newx2;
            break;
        case 10:
            ZerPol=(7.11*newx1*newx1-5.53*newx2*newx2-1.91)*newx1;
            break;
        case 11:
            ZerPol=7.8*rsqr*rsqr-9.36*rsqr+1.9;
            break;
        case 12:
            ZerPol=14.9*(pow(newx1, 4)-pow(newx2, 4))-10*(newx1*newx1-newx2*newx2);
            break;
        case 13:
            ZerPol=(13.2*rsqr-12.4)*newx1*newx2;
            break;
        case 14:
            ZerPol=13*pow(newx1, 4)-2.16*4.39*newx1*newx1*newx2*newx2+12.6*pow(newx2, 4)-1.44*4.39*newx1*newx1-5.97*newx2*newx2+0.71;
            break;
        case 15:
            ZerPol=13.1*(newx1*newx1-newx2*newx2)*newx1*newx2;
            break;
        case 16:
            ZerPol=(18.1*pow(newx1, 4)-36.2*newx1*newx1*newx2*newx2+3.62*pow(newx2, 4)+4.81*newx1*newx1+6.23*newx2*newx2-2.92)*newx2;
            break;
        case 17:
            ZerPol=(11.1*pow(newx1, 4)-22.2*newx1*newx1*newx2*newx2-33.3*pow(newx2, 4)-3.84*newx1*newx1+33*newx2*newx2-2.9)*newx1;            
            break;
        case 18:
            ZerPol=(32.1*pow(newx1, 4)+16.8*newx1*newx1*newx2*newx2+22.7*pow(newx2, 4)-29.2*newx1*newx1-24.6*newx2*newx2+7.21)*newx2;            
            break;
        case 19:
            ZerPol=(34.8*pow(newx1, 4)+16.4*newx1*newx1*newx2*newx2-8.4*pow(newx2, 4)-3.72*newx1*newx1-6.77*newx2*newx2+7.9)*newx1;
            break;
        case 20:
            ZerPol=(14.2*pow(newx1, 4)+3.75*newx1*newx1*newx2*newx2-42.9*pow(newx2, 4)-11.6*newx1*newx1+37*newx2*newx2-5.22)*newx2;            
            break;
        case 21:
            ZerPol=(30*pow(newx1, 4)-21.4*newx1*newx1*newx2*newx2+24.8*pow(newx2, 4)-20.6*newx1*newx1-6.59*newx2*newx2+3.07)*newx1;
            break;
        default:
            throw std::runtime_error("Unidentified Zernike index");
    }

    
    return ZerPol;   
    // return Zernike::ZernikePolynomial(n,m,r,theta)/norm;
}

template<size_t Order>
struct ZernikeAberrations
{
    static std::array<double, (Order + 1) * (Order + 1)> gen_ampls(double sigma, size_t beam_ind);
    static std::array<double, (Order + 1) * (Order + 1)> gen_ampls(std::vector<double> &l_ampls, size_t beam_ind);

    const std::array<double, (Order + 1) * (Order + 1)> m_ampls;
    const double m_beam_radius;
    

    ZernikeAberrations(double sigma, double beam_radius, size_t beam_ind)
        : m_ampls{gen_ampls(sigma, beam_ind)}, m_beam_radius{beam_radius}
    {}

    ZernikeAberrations(std::vector<double> &l_ampls, double beam_radius, size_t beam_ind)
        : m_ampls{gen_ampls(l_ampls, beam_ind)}, m_beam_radius{beam_radius}
    {}

    double operator()(double x1, double x2) const;
    // void load_ampls(std::vector<double> &l_ampls);
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

    double summ=0.0;
    for (size_t i=0; i<m_aberr.size(); i++){
        summ+=m_aberr[i];
    }
    return summ;
    // return std::accumulate(m_aberr.begin(), m_aberr.end(), 0);
}

// template<size_t Order>
// void ZernikeAberrations<Order>::load_ampls(std::vector<double> &l_ampls)
// {
//     size_t l_ampl_ind = 0;
//     for (size_t i = Order; i > 0; i--) {
//         for (size_t j = 0; j < i + 1; j++) {
//             m_ampls[i * Order + j]=l_ampls[l_ampl_ind];
//             l_ampl_ind++;
//         }
//     }
// }

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
    // double ampl_in_mode = dist(rng);
    // double ampl_in_mode = 1;
    for (size_t i = Order; i > 0; i--) {
        for (size_t j = 0; j < i + 1; j++) {
            double ampl_in_mode = dist(rng);
            ampls_sum+=ampl_in_mode;
            m_ampls[i * Order + j] = ampl_in_mode;
            outr /*<< "mode " << i << " " << j << " with coef "*/
                 << ampl_in_mode << '\n';
            mode_num++;
        }
    }
    // double avg=ampls_sum/mode_num;
    // double diffs_sum=0.0;
    // for (size_t i = Order; i > 0; i--) {
    //     for (size_t j = 0; j < i + 1; j++) {
    //         m_diffs.push_back(sqr(m_ampls[i * Order + j]-avg));
    //         diffs_sum+=sqr(m_ampls[i * Order + j]-avg);
    //     }
    // }
    // double sko=sqrt(diffs_sum/mode_num);
    // outr << sko/M_PI*180 << '\n';
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

template<size_t Order>
std::array<double, (Order + 1) * (Order + 1)> ZernikeAberrations<Order>::gen_ampls(std::vector<double> &l_ampls, size_t beam_ind)
{
    std::string fname="poly_ampl_repeat"+std::to_string(beam_ind)+".txt";
    std::ofstream outr(fname);
    std::array<double, (Order + 1) * (Order + 1)> m_ampls;
    size_t l_ampl_ind = 0;
    for (size_t i = Order; i > 0; i--) {
        for (size_t j = 0; j < i + 1; j++) {
            m_ampls[i * Order + j]=l_ampls[l_ampl_ind];
            outr /*<< "mode " << i << " " << j << " with coef "*/
                 << m_ampls[i * Order + j] << '\n';
            l_ampl_ind++;
        }
    }

    return m_ampls;
}



template<size_t MaxInd>
struct ZernikeAberrationsOnSquare
{
    static std::array<double, MaxInd> gen_ampls(double sigma, size_t beam_ind);
    static std::array<double, MaxInd> gen_ampls(std::vector<double> &l_ampls, size_t beam_ind);

    const std::array<double, MaxInd> m_ampls;
    const double m_beam_radius;
    

    ZernikeAberrationsOnSquare(double sigma, double beam_radius, size_t beam_ind)
        : m_ampls{gen_ampls(sigma, beam_ind)}, m_beam_radius{beam_radius}
    {}

    ZernikeAberrationsOnSquare(std::vector<double> &l_ampls, double beam_radius, size_t beam_ind)
        : m_ampls{gen_ampls(l_ampls, beam_ind)}, m_beam_radius{beam_radius}
    {}

    double operator()(double x1, double x2) const;
    // void load_ampls(std::vector<double> &l_ampls);
};

template<size_t MaxInd>
double ZernikeAberrationsOnSquare<MaxInd>::operator()(double x1, double x2) const
{
    std::vector<double> m_aberr;
    for (long j=2; j<MaxInd+1; j++) {
        m_aberr.push_back(ZernikeCartesianOnSquare1(j, x1, x2, m_beam_radius));
        // m_aberr.push_back(ZernikeCartesianOnSquare2(j, x1, x2, m_beam_radius));
        double ampl_in_mode = m_ampls[j];
        m_aberr[j] *= ampl_in_mode;
    }

    double summ=0.0;
    for (size_t i=0; i<m_aberr.size(); i++){
        summ+=m_aberr[i];
    }
    return summ;
    // return std::accumulate(m_aberr.begin(), m_aberr.end(), 0);
}

// template<size_t Order>
// void ZernikeAberrations<Order>::load_ampls(std::vector<double> &l_ampls)
// {
//     size_t l_ampl_ind = 0;
//     for (size_t i = Order; i > 0; i--) {
//         for (size_t j = 0; j < i + 1; j++) {
//             m_ampls[i * Order + j]=l_ampls[l_ampl_ind];
//             l_ampl_ind++;
//         }
//     }
// }

template<size_t MaxInd>
std::array<double, MaxInd> ZernikeAberrationsOnSquare<MaxInd>::gen_ampls(double sigma, size_t beam_ind)
{
    std::string fname="poly_ampl"+std::to_string(beam_ind)+".txt";
    std::ofstream outr(fname);
    std::random_device dev;
    std::mt19937 rng(dev());
    std::array<double, MaxInd> m_ampls;
    std::normal_distribution<double> dist{0, sigma};
    size_t mode_num=0;
    // double ampls_sum=0.0;
    std::vector<double> m_diffs;
    // double ampl_in_mode = dist(rng);
    // double ampl_in_mode = 1;
    for (size_t j = 2; j < MaxInd+1; j++) {
        double ampl_in_mode = dist(rng);
        // ampls_sum+=ampl_in_mode;
        m_ampls[j] = ampl_in_mode;
        outr /*<< "mode " << i << " " << j << " with coef "*/
            << ampl_in_mode << '\n';
    }
      // double avg=ampls_sum/mode_num;
    // double diffs_sum=0.0;
    // for (size_t i = Order; i > 0; i--) {
    //     for (size_t j = 0; j < i + 1; j++) {
    //         m_diffs.push_back(sqr(m_ampls[i * Order + j]-avg));
    //         diffs_sum+=sqr(m_ampls[i * Order + j]-avg);
    //     }
    // }
    // double sko=sqrt(diffs_sum/mode_num);
    // outr << sko/M_PI*180 << '\n';
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

template<size_t MaxInd>
std::array<double, MaxInd> ZernikeAberrationsOnSquare<MaxInd>::gen_ampls(std::vector<double> &l_ampls, size_t beam_ind)
{
    std::string fname="poly_ampl_repeat"+std::to_string(beam_ind)+".txt";
    std::ofstream outr(fname);
    std::array<double, MaxInd> m_ampls;
    size_t l_ampl_ind = 0;

    for (size_t j=2; j < MaxInd+1; j++) {
        m_ampls[j]=l_ampls[l_ampl_ind];
        outr /*<< "mode " << i << " " << j << " with coef "*/
            << m_ampls[j] << '\n';
        l_ampl_ind++;
    }


    return m_ampls;
}


#endif // DISTORTEDBEAM_HPP