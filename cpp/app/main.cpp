#include "stratton-chu/plane-surface.hpp"
#include "stratton-chu/parabolic-surface.hpp"
#include "stratton-chu/elliptic-surface.hpp"
#include "stratton-chu/parallelepiped-volume.hpp"
#include "stratton-chu/distorted-surface.hpp"
#include "stratton-chu/parallel-beam.hpp"
#include "stratton-chu/stratton-chu-field.hpp"
#include "stratton-chu/utils.hpp"
#include "stratton-chu/spec-inverse-fourier.hpp"
#include "stratton-chu/fields-in-region.hpp"
#include "stratton-chu/fields-in-3dregion.hpp"
#include "stratton-chu/distorted-beam.hpp"

#include "stratton-chu/vtk-field-plotter.hpp"
#include "stratton-chu/csv-saver.hpp"
#include "stratton-chu/vtk-saver.hpp"

#include <iostream>
#include "threading.hpp"

#include <cmath>
// #include <fftw3.h>

#include <fstream>
#include <string>

#include <vector>

#include "constants.hpp"
#include "cache.hpp"
#include "static_local_tracker.hpp"

#include <chrono>
#include <random>

#include <zernike_bits/zernike_radial_poly.h>
#include <zernike_bits/zernike_poly.h>

// #include "stdafx.h"
// #include "interpolation.h"

#define WITH_NAME(var) var, #var
#define GET_VARIABLE_NAME(var) (#var)

// const double lambda = 0.000091; // cm = 910 nm
// const double beam_radius = 10; // cm
// const double beam_width = 2 * beam_radius; // cm

//ffff

int main(int argc, const char* argv[])
{
    //Код для монохроматических пучков
    
    //const double F = 20.0; // cm
    // const double phi = M_PI / 3;

    //----------------------------------------------------------------elliptic_mirror_config--------------------------------------------------
    // const double F = 3*beam_width; // cm
    // const double phi = M_PI / 4;


    // double alpha = 40.0/180.0*M_PI;

    // double p = get_p_by_beam_parameters_alpha(alpha, F); // impact parameter по нижней границе
    // double h = p + beam_radius; 
    //-----------------------------------------------------------------------------------------------------------------------------------------

    // std::cout << "Running stratton-chu computation " << argv[1] << ' ' << argv[2] /*<< ' ' << argv[3]<< ' ' << argv[4]*/<<std::endl;
    // int steps_count = 50;

//    double ampl = 0.1 * lambda;
//    double k = 2 * M_PI / (beam_radius / 10);
//    double k = 2 * M_PI / (100 * lambda);

    // double alpha_min = 0.0;
    // double alpha_max = M_PI - phi;
    size_t beams_number = std::atoll(argv[1]);
    // size_t beams_number_con = beams_number*2;
    // //double belt
    // double alpha = M_PI/2;
    // // double p = beam_radius*sqrt(3);
    // double p = beam_radius/tan(M_PI/beams_number);    
    // double F = get_F_by_beam_parameters_alpha_and_p(alpha, p);
    // double h = p + beam_radius; 

    // //double belt cfg 2

    // double alpha = 35.6*M_PI/180;
    // double phi = 43*M_PI/180; 
    // double gamma_med = M_PI/6; 

    double betta, gamma_med, delta, alpha, phic, phi, deltab;
    deltab=0.1;
    // deltab=0.0;
    betta= 2*M_PI/beams_number-deltab;
    // betta= 2*M_PI/beams_number_con-deltab;
    // delta = 0.075;
    delta = 0.1;
    // delta=0.0;
    double forbetta=4;
    // double vert_move=tan(delta)*beam_radius/tan(betta/2);
    double cosalpha, tanphi;
    if (beams_number==4){
        betta= 2*M_PI/beams_number-0.06195919*forbetta;    
        std::cout << tan(betta/2);
        std::cout << std::endl;    
        cosalpha = 0.566051;
        tanphi = 3.08605;
        alpha = acos(cosalpha)+delta;
        // std::cout << alpha;
        // std::cout << std::endl;  
        phi = atan(tanphi)-delta;
        gamma_med=acos(cosalpha);
        std::cout << "Counting 4";
        std::cout << std::endl;

    } else if (beams_number%2==0){        
        gamma_med = asin(tan(betta/2));        
        phic = atan(2*tan(gamma_med)); 
        
        // delta = gamma_med*0.185;
        alpha = gamma_med + delta;
        //alpha = gamma_med;
        phi = phic-delta;
        // phi = phic; 
    // } else if {
    } else {
        gamma_med = asin(sin(betta)/(cos(betta)+cos(betta/2)));      
        phic = atan(2*tan(betta/2)/cos(gamma_med));  
        // delta = 0.1;
        // delta = gamma_med*0.185;
        alpha = gamma_med + delta;
        //alpha = gamma_med;
        phi = phic-delta;
        // phi = phic;      
    }

    // const double F = 3*beam_width; // cm
    const double F = 50.5; // cm
    gamma_med = 9.7323449612/180.0*M_PI; 
    // delta=0.275;
    alpha = 29.7323440615/180.0*M_PI;

    double addlength=17.5307959233;
    double addheight=20.8923888492;





    // double addheight=27.19;
    // double addlength=4.68;

    // double addheight=19.07;
    // double addlength=19.68;

    // double addheight=0.0;
    // double addlength=0.0;

    // double addheight=15.287;
    // double addlength=6.23;

    // double inc_angle = 14.4418/180.0*M_PI;

    // double height_D = sqr(2*30*cos(inc_angle))-4*(30*30-10*10);
    // double height_x = 0.0, height_x1, height_x2;
    
    // if (height_D<0.0) {std::cout << "height check failed as D<0 " <<  std::endl;}
    // else {
    //         height_x1 = (2*30*cos(inc_angle) + sqrt(height_D)) / 2.0;
    //         height_x2 = (2*30*cos(inc_angle) - sqrt(height_D)) / 2.0;
    //         height_x = std::max(height_x1, height_x2);
    //         if (height_x <= 0.0) {std::cout << "height check failed as all x are <0; height doesn't exist " <<  std::endl;}
    // }
    // std::cout << "addheight = " << height_x << std::endl;

    // double add_tan=addlength/addheight;
    // addheight=sqrt(height_x*height_x/(1+add_tan*add_tan));
    // std::cout << "addheight = " << addheight << std::endl;
    // addlength=addheight*add_tan;
    // std::cout << "addlength = " << addlength << std::endl;

    // double addheight=20.89;
    // double addheight=height_x;
    // double addlength=0.0;

      
    // double F = get_F_by_beam_parameters_alpha(alpha, phi, beam_width);
    double p = get_p_by_beam_parameters_alpha(alpha, F); // impact parameter по нижней границе
    double h = p + beam_radius; 

    // // one beam


    // double F = beam_width;
    // double p = get_p_by_beam_parameters_alpha(alpha, F); // impact parameter по нижней границе
    // double h = p + beam_radius; 

    // double x=2*p*tan(alpha-M_PI/6)/(1+sqrt(3)*tan(alpha-M_PI/6));
    // std::cout << x << std::endl; 


    // for(int i = 40; i < steps_count; i++)
    // {
    //     double alpha = alpha_min + (alpha_max - alpha_min) / (steps_count-1) * i;

    // }
    
    
    // std::vector<StrattonChuReflection> refs1;
    // refs1.reserve(beams_number);
    // std::vector<StrattonChuReflection> refs2;
    // refs2.reserve(beams_number);

    std::cout << "Plotting for alpha = " << alpha/M_PI*180; 
    std::cout << " phi = " << phi/M_PI*180; 

    std::cout << "\t Focal length = " << F;
    std::cout << "\t  " << p;
    std::cout << std::endl;

    double sigma = beam_radius;
    double gauss_A = 1;
    // double mirror_radius = 3 * sigma;
    // double mirror_radius = 1.5 * sigma;
    double mirror_radius = 1.1 * sigma;


    // std::vector<ParabolicSurface> mirrors;
    // mirrors.reserve(beams_number);
    std::vector<ParabolicSurface> mirrors1;
    mirrors1.reserve(beams_number);
    std::vector<ParabolicSurface> mirrors2;
    mirrors2.reserve(beams_number);

    std::vector<PlaneSurface> beam_profs1;
    beam_profs1.reserve(beams_number);
    std::vector<PlaneSurface> beam_profs2;
    beam_profs2.reserve(beams_number);


    // std::vector<ParallelBeamAlpha> beams;
    // beams.reserve(beams_number);
    std::vector<ParallelBeamAlpha> beams1;
    beams1.reserve(beams_number);
    std::vector<ParallelBeamAlpha> beams2;
    beams2.reserve(beams_number);

    Position p_focus = {0.0, 0.0, 0.0};
  
    SurfaceRegion region1;
    // region1.x_min = 0.0;
    region1.x_min = h - mirror_radius;
    region1.x_max = h + mirror_radius;
    region1.y_min = 0.0 - mirror_radius;
    region1.y_max = 0.0 + mirror_radius;

    SurfaceRegion region2; //           Почему не минус h ???
    region2.x_min = h - mirror_radius;
    region2.x_max = h + mirror_radius;
    region2.y_min = 0.0 - mirror_radius;
    region2.y_max = 0.0 + mirror_radius;

    // SurfaceRegion region1;
    // region1.x_min = h - beam_radius;
    // region1.x_max = h + beam_radius;
    // region1.y_min = 0.0 - beam_radius;
    // region1.y_max = 0.0 + beam_radius;

    // SurfaceRegion region2; //           Почему не минус h ???
    // region2.x_min = h - beam_radius;
    // region2.x_max = h + beam_radius;
    // region2.y_min = 0.0 - beam_radius;
    // region2.y_max = 0.0 + beam_radius;

    // SurfaceRegion region_profile;
    // region_profile.x_min = h - mirror_radius;
    // region_profile.x_max = h + mirror_radius;
    // region_profile.y_min = 0.0 - mirror_radius;
    // region_profile.y_max = 0.0 + mirror_radius;

    // SurfaceRegion region_profile;
    // region_profile.x_min = h - beam_radius;
    // region_profile.x_max = h + beam_radius;
    // region_profile.y_min = 0.0 - beam_radius;
    // region_profile.y_max = 0.0 + beam_radius;

    SurfaceRegion region_profile;
    region_profile.x_min = 0.0 - beam_radius;
    region_profile.x_max = 0.0 + beam_radius;
    region_profile.y_min = 0.0 - beam_radius;
    region_profile.y_max = 0.0 + beam_radius;

//     // //-----------------------------------------------------------------mirr moves--------------------------------------------------------------------------
//     std::vector<double> mirr_adj;
//     mirr_adj.reserve(2*beams_number);
//     std::random_device dev0;
//     std::mt19937 rng0(dev0());
//     double sigma_mirr=std::atof(argv[2]);

//     std::normal_distribution<double> dist0{0, sigma_mirr};
//     std::string moves_file = "vert_moves_for_"+std::to_string(2*beams_number)+"_beams.txt";
//     // std::string moves_file = "horiz_moves_for_"+std::to_string(2*beams_number)+"_beams.txt";
//     // std::string moves_file = "angz_moves_for_"+std::to_string(2*beams_number)+"_beams.txt";
//     // std::string moves_file = "angy_moves_for_"+std::to_string(2*beams_number)+"_beams.txt";
//     bool rread=false;

//     std::ifstream inm(moves_file);
//     std::string line, firstline;
//     size_t ii=0;   
//     if (inm.is_open()){
//         rread=true;
//         getline (inm,firstline);
//         while (getline (inm,line)) {
//             mirr_adj[ii]=std::stod(line);
//             ii++;
//         }
//         inm.close();
//     } else {
//         std::cout << "Unable to open file" << '\n';
//     }


//     std::ofstream outm(moves_file);
//     if (rread)
//     {
//         std::cout << "meow1" << std::endl;
//         outm << firstline << '\n';
//         for (size_t i=0;i<2*beams_number;i++){
//             outm << mirr_adj[i] << '\n';
//         } 
//     } else {
//         std::cout << "meow2" << std::endl;
//         outm << "for sigma = " <<  sigma_mirr << '\n';
//         for (size_t i=0;i<2*beams_number;i++){
//             mirr_adj[i]=dist0(rng0);
//             // mirr_adj[i]=-beam_radius;
//             outm << mirr_adj[i] << '\n';
//         }
//     }
//     outm.close();

//     double theta0;
// //---------------------------------------------------------------------------end of mirr_moves---------------------------------------------------------------

// // -----------------------------------------------------------------Zernike aberrs--------------------------------------------------------------------------
    std::vector<double> sigma_poly(11);
    sigma_poly[0]=0.1;
    sigma_poly[1]=0.2;
    sigma_poly[2]=0.3;
    sigma_poly[3]=0.4;
    sigma_poly[5]=M_PI/6;
    sigma_poly[6]=M_PI/4;
    sigma_poly[7]=M_PI/3;
    sigma_poly[8]=M_PI/2;
    sigma_poly[9]=M_PI;
    sigma_poly[10]=2*M_PI;
    // sigma_poly[0]=0.075;
    // sigma_poly[1]=0.125;
    
    // sigma_poly[3]=1.0;
    int sig_ind=std::atoi(argv[2]);
    // int sig_ind=3;
    std::cout << sigma_poly[sig_ind]/M_PI*180 << std::endl;

// Zernike standalone aberrs


    std::vector<double> beam_aberr;
    beam_aberr.reserve(2*beams_number);
    std::random_device dev0;
    std::mt19937 rng0(dev0());

    // int zern_ind_n=std::atoi(argv[3]);
    // int zern_ind_m=std::atoi(argv[4]);

    int zern_ind_j=std::atoi(argv[3]);
    // int zern_ind_n=4;
    // int zern_ind_m=2;


    std::normal_distribution<double> dist0{0, sigma_poly[sig_ind]};
    std::string fname="poly_ampl_"+std::to_string(zern_ind_j)+".txt";
    std::ofstream outm(fname);
    outm << "for sigma = " <<  sigma_poly[sig_ind] << '\n';
    for (size_t i=0;i<2*beams_number;i++){
        beam_aberr[i]=dist0(rng0);
        outm << beam_aberr[i] << '\n';
    }
    outm.close();
    
//     // double power_perc_to_side_modes=5;
// // ---------------------------------------------------------------------------end of Zernike aberrs---------------------------------------------------------------


    
    // beams_number=beams_number_con;

    double add_angle;
    if (beams_number%2==1){
        add_angle=M_PI/beams_number;
        // add_angle=0;
    }
    else {
        // add_angle=M_PI/beams_number;
        add_angle=0;

        // // //for con task
        // add_angle=2*M_PI/beams_number;
    }
    
    Vector ax1m1, ax2m1, axFm1;
    Vector ax1m2, ax2m2, axFm2;
    // size_t i_con=0;


// //---------------------------------------------------------------------------------------test PROFILE SEGMENT--------------------------------------
//     for (int zern_ind_j=16; zern_ind_j<22; zern_ind_j++) {
//     // int zern_ind_j=2;

//         ParabolicSurface mirror0(
//             {0.0, 0.0, -F},
//             {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
//             4.0*F, 4.0*F
//         );
//             ParallelBeamAlpha bbeam(lambda, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0},
//                                     //  [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 - h) + sqr(x2)) / (2 * sqr(sigma)) ); },
//                                     // [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
//                                     // [h, gauss_A, sigma, zern=ZernikeAberrations<Order>(sigma_poly[sig_ind], beam_radius*sqrt(2), i)](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*zern(x1-h, x2));},
//                                     //[h, gauss_A, sigma, iw, &beam_aberr, zern_ind_n, zern_ind_m](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[iw]*ZernikeCartesian(zern_ind_n,zern_ind_m,x1-h,x2, beam_radius*sqrt(2)));},
//                                     //  [h, gauss_A, sigma, iw, &beam_aberr, zern_ind_j](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[iw]*ZernikeCartesianOnSquare1(zern_ind_j,x1-h,x2, beam_radius*sqrt(2)));},
//                                     // [h, gauss_A, sigma, &beam_aberr, zern_ind_n, zern_ind_m](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*0.1*ZernikeCartesian(zern_ind_n,zern_ind_m,x1-h,x2, beam_radius*sqrt(2)));},
//                                     // [h, gauss_A, sigma, &beam_aberr, zern_ind_j](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*0.1*ZernikeCartesianOnSquare2(zern_ind_j,x1-h,x2, beam_radius*sqrt(2)));},
//                                     // [gauss_A, zern_ind_n, zern_ind_m](double x1, double x2) { return gauss_A * ZernikeCartesian(zern_ind_n, zern_ind_m, x1, x2, beam_radius*sqrt(2));},
//                                     [gauss_A, zern_ind_j](double x1, double x2) { return gauss_A * ZernikeCartesianOnSquare2(zern_ind_j, x1, x2, beam_radius);},
//                                     // [h](double x1, double x2) { return sqr(x1 - h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
//                                     // [h](double x1, double x2) { return smoothed(sqrt(sqr(x1 - h) + sqr(x2)), beam_radius); },
//                                     [](double, double) { return 0.0; });

//             PlaneSurface beam_prof0({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0});
//             std::string filename_suffixx = "ye_full_"+std::to_string(zern_ind_j);
//             // //plotting mirrors
//             // plot_field_on_given_surface(mirror0, bbeam, region1, 100, "E", "mirror1", filename_suffixx);
//             // std::cout << "mirror0 plotted" << std::endl;

//             //plotting beams
//             plot_field_on_given_surface(beam_prof0, bbeam, region_profile, 300, "E", "beam_profile1", filename_suffixx);
//             std::cout << "beam1 plotted" << std::endl;
//     }
// // //-------------------------------------------------------------------------------------END OF PROFILE SEGMENT------------------------------------


    for (size_t i=0;i<beams_number;i++){
    // for (size_t i=0;i<beams_number;i=i+2){
    // size_t i=0;


//         // //-----------------------------------------double belt, angles_z------------------------------------------------------------------
//         // axFm1={F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)};
//         // ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i), 1.0*cos(gamma_med)};
//         // ax2m1=(ax1m1 % axFm1);
//         // ax2m1 /= ax2m1.norm();
//         // ParabolicSurface mirrorold1(
//         // {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)},
//         // ax1m1,ax2m1, 4.0*F, 4.0*F
//         // );   

//         // theta0=gamma_med+10.0/180.0*M_PI;
//         // // theta0=gamma_med+mirr_adj[i];
//         // axFm1={F*cos(theta0)*cos(2*M_PI/beams_number*i), F*cos(theta0)*sin(2*M_PI/beams_number*i), -F*sin(theta0)};
//         // ax1m1={1.0*sin(theta0)*cos(2*M_PI/beams_number*i), 1.0*sin(theta0)*sin(2*M_PI/beams_number*i), 1.0*cos(theta0)};
//         // ax2m1=(ax1m1 % axFm1);
//         // ax2m1 /= ax2m1.norm();
//         // ParabolicSurface mirrortilted1(
//         // {F*cos(theta0)*cos(2*M_PI/beams_number*i), F*cos(theta0)*sin(2*M_PI/beams_number*i), -F*sin(theta0)},
//         // ax1m1,ax2m1, 4.0*F, 4.0*F
//         // );

//         // Position tilt_r01={F*cos(theta0)*cos(2*M_PI/beams_number*i), F*cos(theta0)*sin(2*M_PI/beams_number*i), -F*sin(theta0)};
//         // Position new_r01=tilt_r01+mirrorold1.point({p,0})-mirrortilted1.point({p,0});

//         // mirrors1.emplace_back(ParabolicSurface(
//         //     {F*cos(theta0)*cos(2*M_PI/beams_number*i), F*cos(theta0)*sin(2*M_PI/beams_number*i), -F*sin(theta0)},
//         //     /*new_r01,*/
//         //     ax1m1,ax2m1, 4.0*F, 4.0*F
//         // ));


//         // axFm2={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
//         // ax1m2={-1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), -1.0*cos(gamma_med)};
//         // ax2m2=(ax1m2 % axFm2);
//         // ax2m2 /= ax2m2.norm();
//         // ParabolicSurface mirrorold2(
//         // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
//         // ax1m2, ax2m2, 4.0*F, 4.0*F
//         // );

//         // theta0=gamma_med+10.0/180.0*M_PI;
//         // // theta0=gamma_med+mirr_adj[i+beams_number];
//         // axFm2={-F*cos(theta0)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(theta0)*sin(2*M_PI/beams_number*i+add_angle), F*sin(theta0)};
//         // ax1m2={-1.0*sin(theta0)*cos(2*M_PI/beams_number*i+add_angle), -1.0*sin(theta0)*sin(2*M_PI/beams_number*i+add_angle), -1.0*cos(theta0)};
//         // ax2m2=(ax1m2 % axFm2);
//         // ax2m2 /= ax2m2.norm();
//         // ParabolicSurface mirrortilted2(
//         // {-F*cos(theta0)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(theta0)*sin(2*M_PI/beams_number*i+add_angle), F*sin(theta0)},
//         // ax1m2, ax2m2, 4.0*F, 4.0*F
//         // );

//         // Position tilt_r02={-F*cos(theta0)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(theta0)*sin(2*M_PI/beams_number*i+add_angle), F*sin(theta0)};
//         // Position new_r02=tilt_r02+mirrorold2.point({p,0})-mirrortilted2.point({p,0});

//         // mirrors2.emplace_back(ParabolicSurface(
//         //     {-F*cos(theta0)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(theta0)*sin(2*M_PI/beams_number*i+add_angle), F*sin(theta0)},
//         //     /*new_r02,*/
//         //     ax1m2,ax2m2, 4.0*F, 4.0*F
//         // ));
        
//         // // axFm1={F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)};
//         // // ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i), 1.0*cos(gamma_med)};
//         // // ax2m1=(ax1m1 % axFm1);
//         // // ax2m1 /= ax2m1.norm();
//         // // axFm2={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
//         // // ax1m2={-1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), -1.0*cos(gamma_med)};
//         // // ax2m2=(ax1m2 % axFm2);
//         // // ax2m2 /= ax2m2.norm();
//         // // //-------------------------------------------------------------------------------------------------------------------------------------------



//         // //-----------------------------------------------double belt, angles_y-----------------------------------------------------------------
//         // double ax1m1_x=1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i);
//         // double ax1m1_y=1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i);
//         // double ax1m1_z=1.0*cos(gamma_med);
//         // axFm1={F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)};
//         // ax1m1={ax1m1_x, ax1m1_y, ax1m1_z};
//         // ax2m1=(ax1m1 % axFm1);
//         // ax2m1 /= ax2m1.norm();
//         // ParabolicSurface mirrorold1(
//         // {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)},
//         // ax1m1,ax2m1, 4.0*F, 4.0*F
//         // );   

//         // // Position fornewax1=mirrorold1.point({p+beam_radius,0})-mirrorold1.point({p,0});
//         // // Vector newax1={fornewax1[0],fornewax1[1],fornewax1[2]};
//         // // newax1 /= newax1.norm();
//         // theta0=mirr_adj[i];
//         // Vector axFm1_new;
//         // double M1[3][3] = {
//         //     {cos(theta0)+(1-cos(theta0))*ax1m1_x*ax1m1_x,(1-cos(theta0))*ax1m1_x*ax1m1_y-sin(theta0)*ax1m1_z,(1-cos(theta0))*ax1m1_x*ax1m1_z+sin(theta0)*ax1m1_y},
//         //     {(1-cos(theta0))*ax1m1_x*ax1m1_y+sin(theta0)*ax1m1_z, cos(theta0)+(1-cos(theta0))*ax1m1_y*ax1m1_y,(1-cos(theta0))*ax1m1_y*ax1m1_z-sin(theta0)*ax1m1_x},
//         //     {(1-cos(theta0))*ax1m1_x*ax1m1_z-sin(theta0)*ax1m1_y,(1-cos(theta0))*ax1m1_y*ax1m1_z+sin(theta0)*ax1m1_x,cos(theta0)+(1-cos(theta0))*ax1m1_z*ax1m1_z}
//         // };

//         // for (size_t j=0; j<3;j++){
//         //     axFm1_new[j]=M1[j][0]*axFm1[0]+M1[j][1]*axFm1[1]+M1[j][2]*axFm1[2];
//         //     // ax1m1_new[j]=M1[j][0]*ax1m1[0]+M1[j][1]*ax1m1[1]+M1[j][2]*ax1m1[2];
//         // }
//         // Vector ax2m1_new;
//         // ax2m1_new=(ax1m1 % axFm1_new);
//         // ax2m1_new /= ax2m1_new.norm();
//         // ParabolicSurface mirrortilted1(
//         //     {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)},
//         //     ax1m1,ax2m1_new, 4.0*F, 4.0*F
//         // );

//         // Position old_r01={F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)};
//         // Position new_r01=old_r01+mirrorold1.point({p,0})-mirrortilted1.point({p,0});

//         // mirrors1.emplace_back(ParabolicSurface(
//         //     new_r01,
//         //     ax1m1,ax2m1_new, 4.0*F, 4.0*F
//         // ));

//         // // mirrors1.emplace_back(ParabolicSurface(
//         // //     {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)},
//         // //     ax1m1,ax2m1_new, 4.0*F, 4.0*F
//         // // ));


//         // double ax1m2_x=-1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle);
//         // double ax1m2_y=-1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle);
//         // double ax1m2_z=-1.0*cos(gamma_med);
//         // axFm2={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
//         // ax1m2={ax1m2_x, ax1m2_y, ax1m2_z};
//         // ax2m2=(ax1m2 % axFm2);
//         // ax2m2 /= ax2m2.norm();
//         // ParabolicSurface mirrorold2(
//         // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
//         // ax1m2, ax2m2, 4.0*F, 4.0*F
//         // );

//         // // Position fornewax2=mirrorold2.point({p+beam_radius,0})-mirrorold2.point({p,0});
//         // // Vector newax2={fornewax2[0],fornewax2[1],fornewax2[2]};
//         // // newax2 /= newax2.norm();
//         // theta0=mirr_adj[i+beams_number];
//         // Vector axFm2_new;
//         // double M2[3][3] = {
//         //     {cos(theta0)+(1-cos(theta0))*ax1m2_x*ax1m2_x,(1-cos(theta0))*ax1m2_x*ax1m2_y-sin(theta0)*ax1m2_z,(1-cos(theta0))*ax1m2_x*ax1m2_z+sin(theta0)*ax1m2_y},
//         //     {(1-cos(theta0))*ax1m2_x*ax1m2_y+sin(theta0)*ax1m2_z, cos(theta0)+(1-cos(theta0))*ax1m2_y*ax1m2_y,(1-cos(theta0))*ax1m2_y*ax1m2_z-sin(theta0)*ax1m2_x},
//         //     {(1-cos(theta0))*ax1m2_x*ax1m2_z-sin(theta0)*ax1m2_y,(1-cos(theta0))*ax1m2_y*ax1m2_z+sin(theta0)*ax1m2_x,cos(theta0)+(1-cos(theta0))*ax1m2_z*ax1m2_z}
//         // };

//         // for (size_t j=0; j<3;j++){
//         //     axFm2_new[j]=M2[j][0]*axFm2[0]+M2[j][1]*axFm2[1]+M2[j][2]*axFm2[2];
//         //     // ax1m2_new[j]=M1[j][0]*ax1m2[0]+M1[j][1]*ax1m2[1]+M1[j][2]*ax1m2[2];
//         // }
//         // Vector ax2m2_new;
//         // ax2m2_new=(ax1m2 % axFm2_new);
//         // ax2m2_new /= ax2m2_new.norm();
//         // ParabolicSurface mirrortilted2(
//         // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
//         // ax1m2, ax2m2_new, 4.0*F, 4.0*F
//         // );

//         // Position old_r02={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
//         // Position new_r02=old_r02+mirrorold2.point({p,0})-mirrortilted2.point({p,0});

//         // mirrors2.emplace_back(ParabolicSurface(
//         //     new_r02,
//         //     ax1m2, ax2m2_new, 4.0*F, 4.0*F
//         // ));

//         // // mirrors2.emplace_back(ParabolicSurface(
//         // //     {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
//         // //     ax1m2, ax2m2_new, 4.0*F, 4.0*F
//         // // ));
//         // //----------------------------------------------------------------------------------------------------------------------


        //-------------------------------------------------------------------double belt, vert+horiz+horiz2-----------------------------------------------------
        axFm1={F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)};
        ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i), 1.0*cos(gamma_med)};
        ax2m1=(ax1m1 % axFm1);
        ax2m1 /= ax2m1.norm();
        mirrors1.emplace_back(ParabolicSurface(
            // {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)},
            // {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)-beam_radius},
            // {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)+mirr_adj[i]},
            // {(F+mirr_adj[i]/cos(gamma_med))*cos(gamma_med)*cos(2*M_PI/beams_number*i), (F+mirr_adj[i]/cos(gamma_med))*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)},
            {(F*cos(gamma_med)+addlength)*cos(2*M_PI/beams_number*i), (F*cos(gamma_med)+addlength)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)+addheight},
            // {F*cos(gamma_med)*cos(2*M_PI/beams_number*i)-sin(2*M_PI/beams_number*i)*mirr_adj[i], F*cos(gamma_med)*sin(2*M_PI/beams_number*i)+cos(2*M_PI/beams_number*i)*mirr_adj[i], -F*sin(gamma_med)},
            ax1m1,ax2m1, 4.0*F, 4.0*F
        ));

        // ParabolicSurface mirror1(
        //     {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i),-F*sin(gamma_med)},
        //     ax1m1,
        //     ax2m1,
        //     4.0*F, 4.0*F
        // );


        axFm2={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
        ax1m2={-1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), -1.0*cos(gamma_med)};
        ax2m2=(ax1m2 % axFm2);
        ax2m2 /= ax2m2.norm();
        mirrors2.emplace_back(ParabolicSurface(
            // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
            // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)+beam_radius},
            // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)-mirr_adj[i+beams_number]},
            // {-(F+mirr_adj[i+beams_number]/cos(gamma_med))*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -(F+mirr_adj[i+beams_number]/cos(gamma_med))*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
            {-(F*cos(gamma_med)+addlength)*cos(2*M_PI/beams_number*i+add_angle), -(F*cos(gamma_med)+addlength)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)-addheight},
            // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle)+sin(2*M_PI/beams_number*i+add_angle)*mirr_adj[i+beams_number], -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle)-cos(2*M_PI/beams_number*i+add_angle)*mirr_adj[i+beams_number], F*sin(gamma_med)},
            ax1m2, ax2m2, 4.0*F, 4.0*F
        ));

        // ParabolicSurface mirror2(
        //     {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
        //     ax1m2,
        //     ax2m2,
        //     4.0*F, 4.0*F
        // );
        // //------------------------------------------------------------------------------------------------

        // // ParallelBeamAlpha beam1(lambda, {0.0, 0.0, 0.0}, ax1m1, -ax2m1,
        // //                         [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
        // //                         [](double, double) { return 0.0; });

        

        beams1.emplace_back(ParallelBeamAlpha(lambda, {0.0, 0.0, 0.0}, /*new_r01,*/ ax1m1, -ax2m1,
        // beams1.emplace_back(ParallelBeamAlpha(lambda, {addlength*cos(2*M_PI/beams_number*i), addlength*sin(2*M_PI/beams_number*i), addheight}, /*new_r01,*/ ax1m1, -ax2m1,
                                //  [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 - h) + sqr(x2)) / (2 * sqr(sigma)) ); },
                                [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
                                // [h, gauss_A, sigma, zern=ZernikeAberrations<Order>(sigma_poly[sig_ind], beam_radius*sqrt(2), i)](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*zern(x1-h, x2));},
                                // [h, gauss_A, sigma, i, &beam_aberr, zern_ind_n, zern_ind_m](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i]*ZernikeCartesian(zern_ind_n,zern_ind_m,x1-h,x2, beam_radius*sqrt(2)));},
                                // [h, gauss_A, sigma, zern=ZernikeAberrations<MaxInd>(sigma_poly[sig_ind], beam_radius*sqrt(2), i)](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*zern(x1-h, x2));},
                                // [h, gauss_A, sigma, i, &beam_aberr, zern_ind_j](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i]*ZernikeCartesianOnSquare1(zern_ind_j,x1-h,x2, beam_radius*sqrt(2)));},
                                //[h](double x1, double x2) { return sqr(x1 - h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
                                // [h](double x1, double x2) { return smoothed(sqrt(sqr(x1 - h) + sqr(x2)), beam_radius); },
                                [](double, double) { return 0.0; }));

        ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), 1.0*cos(gamma_med)};
        // // ax1m1={1.0*sin(theta0)*cos(2*M_PI/beams_number*i+add_angle), 1.0*sin(theta0)*sin(2*M_PI/beams_number*i+add_angle), 1.0*cos(theta0)};

        // ParallelBeamAlpha beam2(lambda, {0.0, 0.0, 0.0}, ax1m1, ax2m2,
        //                         [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
        //                         [](double, double) { return 0.0; });

        beams2.emplace_back(ParallelBeamAlpha(lambda, {0.0, 0.0, 0.0},/*new_r02,*/ ax1m1, ax2m2,
        // beams2.emplace_back(ParallelBeamAlpha(lambda, {-addlength*cos(2*M_PI/beams_number*i+add_angle), -addlength*sin(2*M_PI/beams_number*i+add_angle), -addheight},/*new_r02,*/ ax1m1, ax2m2,
                                //  [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 + h) + sqr(x2)) / (2 * sqr(sigma)) ); },
                                [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
                                // [h, gauss_A, sigma, zern=ZernikeAberrations<Order>(sigma_poly[sig_ind], beam_radius*sqrt(2), i+beams_number)](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*zern(x1+h, x2)); },
                                // [h, gauss_A, sigma, i, beams_number, &beam_aberr, zern_ind_n, zern_ind_m](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i+beams_number]*ZernikeCartesian(zern_ind_n,zern_ind_m,x1+h,x2, beam_radius*sqrt(2))); },
                                // [h, gauss_A, sigma, zern=ZernikeAberrations<MaxInd>(sigma_poly[sig_ind], beam_radius*sqrt(2), i+beams_number)](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*zern(x1+h, x2)); },
                                // [h, gauss_A, sigma, i, beams_number, &beam_aberr, zern_ind_j](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i+beams_number]*ZernikeCartesianOnSquare1(zern_ind_j,x1+h,x2, beam_radius*sqrt(2))); },
                                //[h](double x1, double x2) { return sqr(x1 + h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
                                // [h](double x1, double x2) { return smoothed(sqrt(sqr(x1 + h) + sqr(x2)), beam_radius); },
                                [](double, double) { return 0.0; }));
        ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i), 1.0*cos(gamma_med)};


// //---------------------------------------------------------------straight_con---------------------------------------------------------------------------
//         // //straight_con
//         // axFm1={F*cos(2*M_PI/beams_number*i), F*sin(2*M_PI/beams_number*i), 0.0};
//         // ax1m1={0.0, 0.0, 1.0};
//         // ax2m1=(ax1m1 % axFm1);
//         // ax2m1 /= ax2m1.norm();
//         // mirrors1.emplace_back(ParabolicSurface(
//         //     {F*cos(2*M_PI/beams_number*i), F*sin(2*M_PI/beams_number*i), 0.0},
//         //     // {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)-beam_radius},
//         //     // {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)+mirr_adj[i]},
//         //     // {(F+mirr_adj[i]/cos(gamma_med))*cos(gamma_med)*cos(2*M_PI/beams_number*i), (F+mirr_adj[i]/cos(gamma_med))*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)},
//         //     // {F*cos(gamma_med)*cos(2*M_PI/beams_number*i)-sin(2*M_PI/beams_number*i)*mirr_adj[i], F*cos(gamma_med)*sin(2*M_PI/beams_number*i)+cos(2*M_PI/beams_number*i)*mirr_adj[i], -F*sin(gamma_med)},
//         //     ax1m1,ax2m1, 4.0*F, 4.0*F
//         // ));

//         // ParabolicSurface mirror1(
//         //     {F*cos(2*M_PI/beams_number*i), F*sin(2*M_PI/beams_number*i), 0.0},
//         //     ax1m1,ax2m1, 4.0*F, 4.0*F
//         // );


//         // axFm2={-F*cos(2*M_PI/beams_number*i+add_angle), -F*sin(2*M_PI/beams_number*i+add_angle), 0.0};
//         // ax1m2={0.0, 0.0, -1.0};
//         // ax2m2=(ax1m2 % axFm2);
//         // ax2m2 /= ax2m2.norm();
//         // mirrors2.emplace_back(ParabolicSurface(
//         //     {-F*cos(2*M_PI/beams_number*i+add_angle), -F*sin(2*M_PI/beams_number*i+add_angle), 0.0},
//         //     // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)+beam_radius},
//         //     // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)-mirr_adj[i+beams_number]},
//         //     // {-(F+mirr_adj[i+beams_number]/cos(gamma_med))*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -(F+mirr_adj[i+beams_number]/cos(gamma_med))*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
//         //     // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle)+sin(2*M_PI/beams_number*i+add_angle)*mirr_adj[i+beams_number], -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle)-cos(2*M_PI/beams_number*i+add_angle)*mirr_adj[i+beams_number], F*sin(gamma_med)},
//         //     ax1m2, ax2m2, 4.0*F, 4.0*F
//         // ));

//         // ParabolicSurface mirror2(
//         //     {-F*cos(2*M_PI/beams_number*i+add_angle), -F*sin(2*M_PI/beams_number*i+add_angle), 0.0},
//         //     ax1m2, ax2m2, 4.0*F, 4.0*F
//         // );

//         // beams1.emplace_back(ParallelBeamAlpha(lambda, {0.0, 0.0, 0.0}, /*{0.0, 0.0, -beam_radius},*/ ax1m1, -ax2m1,
//         //                         //  [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 - h) + sqr(x2)) / (2 * sqr(sigma)) ); },
//         //                         [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
//         //                         // [h, gauss_A, sigma, zern=ZernikeAberrations<Order>(sigma_poly[sig_ind], beam_radius*sqrt(2), i)](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*zern(x1-h, x2));},
//         //                         // [h, gauss_A, sigma, i, &beam_aberr, zern_ind_n, zern_ind_m](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i]*ZernikeCartesian(zern_ind_n,zern_ind_m,x1-h,x2, beam_radius*sqrt(2)));},
//         //                         // [h, gauss_A, sigma, zern=ZernikeAberrationsOnSquare<MaxInd>(sigma_poly[sig_ind], beam_radius*sqrt(2), i)](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*zern(x1-h, x2));},
//         //                         // [h, gauss_A, sigma, i, &beam_aberr, zern_ind_j](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i]*ZernikeCartesianOnSquare1(zern_ind_j,x1-h,x2, beam_radius*sqrt(2)));},
//         //                         //[h](double x1, double x2) { return sqr(x1 - h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
//         //                         // [h](double x1, double x2) { return smoothed(sqrt(sqr(x1 - h) + sqr(x2)), beam_radius); },
//         //                         [](double, double) { return 0.0; }));

//         // ParallelBeamAlpha beam1(lambda, {0.0, 0.0, 0.0}, /*{0.0, 0.0, -beam_radius},*/ ax1m1, -ax2m1,
//         //                         [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
//         //                         [](double, double) { return 0.0; });

//         // ax1m1={0.0, 0.0, 1.0};

//         // beams2.emplace_back(ParallelBeamAlpha(lambda, {0.0, 0.0, 0.0},/*{0.0, 0.0, beam_radius},*/ ax1m1, ax2m2,
//         //                         //  [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 + h) + sqr(x2)) / (2 * sqr(sigma)) ); },
//         //                         [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
//         //                         // [h, gauss_A, sigma, zern=ZernikeAberrations<Order>(sigma_poly[sig_ind], beam_radius*sqrt(2), i+beams_number)](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*zern(x1+h, x2)); },
//         //                         // [h, gauss_A, sigma, i, beams_number, &beam_aberr, zern_ind_n, zern_ind_m](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i+beams_number]*ZernikeCartesian(zern_ind_n,zern_ind_m,x1+h,x2, beam_radius*sqrt(2))); },
//         //                         // [h, gauss_A, sigma, zern=ZernikeAberrationsOnSquare<MaxInd>(sigma_poly[sig_ind], beam_radius*sqrt(2), i+beams_number)](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*zern(x1+h, x2)); },
//         //                         // [h, gauss_A, sigma, i, beams_number, &beam_aberr, zern_ind_j](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i+beams_number]*ZernikeCartesianOnSquare1(zern_ind_j,x1+h,x2, beam_radius*sqrt(2))); },

//         //                         //[h](double x1, double x2) { return sqr(x1 + h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
//         //                         // [h](double x1, double x2) { return smoothed(sqrt(sqr(x1 + h) + sqr(x2)), beam_radius); },
//         //                         [](double, double) { return 0.0; }));

//         // ParallelBeamAlpha beam2(lambda, {0.0, 0.0, 0.0},/*{0.0, 0.0, beam_radius},*/ ax1m1, ax2m2,
//         //                         [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
//         //                         [](double, double) { return 0.0; });                        
// //-------------------------------------------------------------------------------------------------------------------------------------------------------


        // refs1.emplace_back(mirrors1[i], beams1[i], region1);
        // refs2.emplace_back(mirrors2[i], beams2[i], region2);
        std::string filename_suffix = std::to_string(i);
        // std::string filename_suffix = std::to_string(i_con);

//         beam_profs1.emplace_back(PlaneSurface({0.0, 0.0, 0.0}, ax1m1, -ax2m1));
//         // PlaneSurface beam_prof1({0.0, 0.0, 0.0}, ax1m1, -ax2m1);
//         beam_profs2.emplace_back(PlaneSurface({0.0, 0.0, 0.0}, ax1m2, ax2m2));

        // //plotting mirrors
        // plot_field_on_given_surface(mirror1, beam1, region1, 100, "E", "mirror1", filename_suffix);
        // std::cout << "mirror1 plotted" << std::endl;
        // plot_field_on_given_surface(mirror2, beam2, region2, 100, "E", "mirror2", filename_suffix);
        // std::cout << "mirror2 plotted" << std::endl;
        plot_field_on_given_surface(mirrors1[i], beams1[i], region1, 100, "E", "mirror1", filename_suffix);
        std::cout << "mirror1 plotted" << std::endl;
        plot_field_on_given_surface(mirrors2[i], beams2[i], region2, 100, "E", "mirror2", filename_suffix);
        std::cout << "mirror2 plotted" << std::endl;

//         // //plotting mirrors (skip con)
//         // plot_field_on_given_surface(mirrors1[i_con], beams1[i_con], region1, 100, "E", "mirror1", filename_suffix);
//         // std::cout << "mirror1 plotted" << std::endl;
//         // plot_field_on_given_surface(mirrors2[i_con], beams2[i_con], region2, 100, "E", "mirror2", filename_suffix);
//         // std::cout << "mirror2 plotted" << std::endl;

//         // //plotting beams
//         // plot_field_on_given_surface(beam_prof1, beam1, region_profile, 300, "E", "beam_profile1", filename_suffix);
//         // std::cout << "beam1 plotted" << std::endl;
//         // plot_field_on_given_surface(beam_profs2[i], beams2[i], region_profile, 300, "E", "beam_profile2", filename_suffix);
//         // std::cout << "beam2 plotted" << std::endl;

//         // i_con++;
    }
//     // beams_number = std::atoll(argv[1]);
    std::cout << "mirror plotting ended" << std::endl;


// // //---------------------------------------------------------------------------------------testing elliptic mirror--------------------------------------

// double testella = 30.0, testellb = 25.0;

// Vector testelax1={1.0, 0.0, 0.0};
// Vector testelax2={0.0, 1.0, 0.0};
// double testellc = sqrt(1-sqr(testellb/testella))*testella;

// std::cout << "c " << testellc << std::endl;

// Vector testelax3=testelax1 % testelax2;

// std::cout << "testelax3 " << testelax3[0] << " " << testelax3[1] << " " << testelax3[2] <<  std::endl;

// EllipticSurface testell_mirror(
//     {0.0, 0.0, 0.0},
//     testelax1, testelax2, testella, testellb, testellb
// );

// PlaneSurface testfocus1({-testellc, 0.0, 0.0}, testelax1, testelax2);
// PlaneSurface testfocus2({testellc, 0.0, 0.0}, testelax1, testelax2);

// ParallelBeamAlpha testbeam(lambda, {-testellc, 0.0, 0.0}, testelax1, testelax2,
//                             [gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1,2) + pow(x2,2)) / (2 * pow(1*lambda, 2)) ); },
//                             [](double, double) { return 0.0; });

// SurfaceRegion testell_beam_region;

// testell_beam_region.x_min = -testellc - 6.0*lambda;
// testell_beam_region.x_max = -testellc + 6.0*lambda;
// testell_beam_region.y_min = 0.0 - 6.0*lambda;
// testell_beam_region.y_max = 0.0 + 6.0*lambda;


// // testell_beam_region.x_min = -testellc - 1.0;
// // testell_beam_region.x_max = -testellc + 1.0;
// // testell_beam_region.y_min = 0.0 - 1.0;
// // testell_beam_region.y_max = 0.0 + 1.0;

// Position focus2={testellc, 0.0, 0.0};
// Position testmirr_cent={-testellc, 0.0, testellb*sqrt(1-testellc*testellc/testella/testella)};

// Vector testaxtoF2=focus2-testmirr_cent;
// testaxtoF2 /= testaxtoF2.norm();
// std::cout << "testaxtoF2 " << testaxtoF2[0] << " " << testaxtoF2[1] << " " << testaxtoF2[2] <<  std::endl;

// Vector testaxF21=-(testaxtoF2 % testelax2);
// testaxF21 /= testaxF21.norm();
// std::cout << "focal testax1 " << testaxF21[0] << " " << testaxF21[1] << " " << testaxF21[2] <<  std::endl;

// PlaneSurface focal_zone_norm({testellc, 0.0, 0.0}, testaxF21, testelax2);

// // SurfaceRegion testfocal_region1;
// // testfocal_region1.x_min = 0.0 - 2.0;
// // testfocal_region1.x_max = 0.0 + 2.0;
// // testfocal_region1.y_min = 0.0 - 2.0;
// // testfocal_region1.y_max = 0.0 + 2.0;

// SurfaceRegion testfocal_region1;
// testfocal_region1.x_min = 0.0 - 5.0*lambda;
// testfocal_region1.x_max = 0.0 + 5.0*lambda;
// testfocal_region1.y_min = 0.0 - 5.0*lambda;
// testfocal_region1.y_max = 0.0 + 5.0*lambda;
// plot_field_on_given_surface(testfocus1, testbeam, testfocal_region1, 100, "E", "ELL_mirror_focus1", "stepTT");
// plot_field_on_given_surface_non_threaded(testell_mirror, testbeam, testell_beam_region, 100, "E", "ELL_mirror", "stepTT");

// StrattonChuReflection testrefl(testell_mirror, testbeam, testell_beam_region);
// plot_field_on_given_surface(focal_zone_norm, testrefl, testfocal_region1, 100, "E", "ELL_mirror_focus2", "stepTT");

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
    SurfaceRegion region1g;
    region1g.x_min = h - mirror_radius;
    region1g.x_max = h + mirror_radius;
    region1g.y_min = 0.0 - mirror_radius;
    region1g.y_max = 0.0 + mirror_radius;

    Position main_focus = {0.0, 0.0, 0.0};
    Position second_focus = {0.0, 0.0, 0.0};

    Vector elax3, axPM2;

    Vector elax1, elax1_copy;

    double ell_a=19.9352856928, ell_b=14.5416792797;

    double ell_c, theta_in, mid_D, mid_x, mid_x1, mid_x2, d_alpha;

    Position ell_ax_cent;

    Vector2D MID2d;
    MID2d[0]=h;
    MID2d[1]=0;

    Position MID3d;

    Vector axtoF;
    Position ell_point_glob, ell_point_loc, new_ell_point_loc;

    Vector axF11;

    SurfaceRegion focal_region;
    // focal_region.x_min = 0.0 - 3 * lambda;
    // focal_region.x_max = 0.0 + 3 * lambda;
    // focal_region.y_min = 0.0 - 3 * lambda;
    // focal_region.y_max = 0.0 + 3 * lambda;
    focal_region.x_min = 0.0 - 5 * lambda;
    focal_region.x_max = 0.0 + 5 * lambda;
    focal_region.y_min = 0.0 - 5 * lambda;
    focal_region.y_max = 0.0 + 5 * lambda;

    SurfaceRegion focal_regionELL2;
    focal_regionELL2.x_min = 0.0 - 2 * lambda;
    focal_regionELL2.x_max = 0.0 + 2 * lambda;
    focal_regionELL2.y_min = 0.0 - 2 * lambda;
    focal_regionELL2.y_max = 0.0 + 2 * lambda;

    std::vector<ParabolicSurface> parmirrors1;
    parmirrors1.reserve(beams_number);
    std::vector<ParabolicSurface> parmirrors2;
    parmirrors2.reserve(beams_number);

    std::vector<EllipticSurface> ellmirrors1;
    ellmirrors1.reserve(beams_number);
    std::vector<EllipticSurface> ellmirrors2;
    ellmirrors2.reserve(beams_number);

    std::vector<PlaneSurface> firstfocus1;
    firstfocus1.reserve(beams_number);
    std::vector<PlaneSurface> firstfocus2;
    firstfocus2.reserve(beams_number);

    std::vector<ParallelBeamAlpha> incbeams1;
    incbeams1.reserve(beams_number);
    std::vector<ParallelBeamAlpha> incbeams2;
    incbeams2.reserve(beams_number);

    std::vector<StrattonChuReflection> firstrefs1;
    firstrefs1.reserve(beams_number);
    std::vector<StrattonChuReflection> firstrefs2;
    firstrefs2.reserve(beams_number);

    std::vector<StrattonChuReflection> focus1_pass1;
    focus1_pass1.reserve(beams_number);
    std::vector<StrattonChuReflection> focus1_pass2;
    focus1_pass2.reserve(beams_number);

    std::vector<StrattonChuReflectionLowPrec> secondrefs1;
    secondrefs1.reserve(beams_number);
    std::vector<StrattonChuReflectionLowPrec> secondrefs2;
    secondrefs2.reserve(beams_number);

    std::vector<SurfaceRegion> ellregs1;
    ellregs1.reserve(beams_number);
    std::vector<SurfaceRegion> ellregs2;
    ellregs2.reserve(beams_number);

    Vector2D low2d;
    low2d[0]=p;
    low2d[1]=0;

    Vector2D high2d;
    high2d[0]=p+beam_width;
    high2d[1]=0;

    Position low3d, high3d;
    
    Vector lowtoF, hightoF;

    double lowtomid, new_lowtomid, low_D, low_x, low_x1, low_x2;
   
    double hightomid, hightolow, new_hightomid, high_D, high_x, high_x1, high_x2;
    
    Vector midtoF, Ftoellmid;
    double ell_mirr_length;



// // //-----------------------------------------------------------------------------------------------------------------------------------------------------
    for (size_t i=0;i<beams_number;i++){
    // size_t i=0;

    axFm1={F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)};
    ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i), 1.0*cos(gamma_med)};
    ax2m1=(ax1m1 % axFm1);
    ax2m1 /= ax2m1.norm();
    parmirrors1.emplace_back(ParabolicSurface(
            {(F*cos(gamma_med)+addlength)*cos(2*M_PI/beams_number*i), (F*cos(gamma_med)+addlength)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)+addheight},
            ax1m1,ax2m1, 4.0*F, 4.0*F
    ));
    std::cout << "ax2m1 " << ax2m1[0] << " " << ax2m1[1] << " " << ax2m1[2] <<  std::endl;
    // ParabolicSurface mirror(
    //     {(F*cos(gamma_med)+addlength)*cos(2*M_PI/beams_number*i), (F*cos(gamma_med)+addlength)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)+addheight},
    //     ax1m1,
    //     ax2m1,
    //     4.0*F, 4.0*F
    // );
    incbeams1.emplace_back(ParallelBeamAlpha(lambda, {addlength*cos(2*M_PI/beams_number*i), addlength*sin(2*M_PI/beams_number*i), addheight}, /*new_r01,*/ ax1m1, -ax2m1,
                                //  [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 - h) + sqr(x2)) / (2 * sqr(sigma)) ); },
                                // [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
                                // [h, gauss_A, sigma, zern=ZernikeAberrations<Order>(sigma_poly[sig_ind], beam_radius*sqrt(2), i)](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*zern(x1-h, x2));},
                                // [h, gauss_A, sigma, i, &beam_aberr, zern_ind_n, zern_ind_m](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i]*ZernikeCartesian(zern_ind_n,zern_ind_m,x1-h,x2, beam_radius*sqrt(2)));},
                                [h, gauss_A, sigma, i, &beam_aberr, zern_ind_j](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i]*ZernikeCartesianOnSquare1(zern_ind_j,x1-h,x2, beam_radius));},
                                //[h](double x1, double x2) { return sqr(x1 - h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
                                // [h](double x1, double x2) { return smoothed(sqrt(sqr(x1 - h) + sqr(x2)), beam_radius); },
                                [](double, double) { return 0.0; }));

    // ParallelBeamAlpha beam(lambda, {addlength*cos(2*M_PI/beams_number*i), addlength*sin(2*M_PI/beams_number*i), addheight}, /*new_r01,*/ ax1m1, -ax2m1,                           
    //                       [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
    //                       [](double, double) { return 0.0; });



    // axFm2={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
    // ax1m2={-1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), -1.0*cos(gamma_med)};
    // ax2m2=(ax1m2 % axFm2);
    // ax2m2 /= ax2m2.norm();

    // ParabolicSurface mirror(
    //     {-(F*cos(gamma_med)+addlength)*cos(2*M_PI/beams_number*i+add_angle), -(F*cos(gamma_med)+addlength)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)-addheight},
    //     ax1m2,
    //     ax2m2,
    //     4.0*F, 4.0*F
    // );
    // ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), 1.0*cos(gamma_med)};
    // ParallelBeamAlpha beam(lambda, {-addlength*cos(2*M_PI/beams_number*i+add_angle), -addlength*sin(2*M_PI/beams_number*i+add_angle), -addheight}, /*new_r01,*/ ax1m1, ax2m2,                          
    //                       [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
    //                       [](double, double) { return 0.0; });
    

    
    p_focus = {addlength*cos(2*M_PI/beams_number*i), addlength*sin(2*M_PI/beams_number*i), addheight};
    elax3 = -ax2m1;
    axPM2 = ax2m1;

    // p_focus = {-addlength*cos(2*M_PI/beams_number*i+add_angle), -addlength*sin(2*M_PI/beams_number*i+add_angle), -addheight};
    // Vector elax3 = -ax2m2;
    // Vector axPM2 = ax2m2;


    elax1=main_focus-p_focus;
    elax1_copy=main_focus-p_focus;
    elax1 /= elax1.norm();


    std::cout << "elax1 " << elax1[0] << " " << elax1[1] << " " << elax1[2] <<  std::endl;
    std::cout << "elax3 " << elax3[0] << " " << elax3[1] << " " << elax3[2] <<  std::endl;

    ell_c=elax1_copy.norm()/2;
    ell_ax_cent = p_focus + elax1 * ell_c;
    std::cout << "ell_ax_cent " << ell_ax_cent[0] << " " << ell_ax_cent[1] << " " << ell_ax_cent[2] <<  std::endl;
    

    ellmirrors1.emplace_back(EllipticSurface(ell_ax_cent, elax1, elax3, ell_a, ell_b, ell_b)); 


  
    std::cout << "plotting parab mirror" << std::endl;

    std::string filename_suffixST1 = "step1_"+std::to_string(i);

    plot_field_on_given_surface(parmirrors1[i], incbeams1[i], region1g, 111, "E", "mirror_par1", filename_suffixST1);

    std::cout << "parab mirror plotted" << std::endl;
    
    firstrefs1.emplace_back(parmirrors1[i], incbeams1[i], region1g);

    
    MID3d=parmirrors1[i].point(MID2d);
    std::cout << "mirror mid " << MID3d[0] << " " << MID3d[1] << " " << MID3d[2] <<  std::endl;

    std::cout << "p_focus " << p_focus[0] << " " << p_focus[1] << " " << p_focus[2] <<  std::endl;

    
    axtoF=-(MID3d-p_focus);
    axtoF /= axtoF.norm();
    std::cout << "axtoF " << axtoF[0] << " " << axtoF[1] << " " << axtoF[2] <<  std::endl;

    theta_in = acos(elax1*axtoF);

    std::cout << "theta_in "<< theta_in/M_PI*180 <<  std::endl;

    mid_D = sqr(ell_b*ell_a)*(sqr(tan(theta_in))*(ell_a*ell_a-ell_c*ell_c)+ell_b*ell_b);
    mid_x = 0.0;
    
    if (mid_D<0.0) {std::cout << "mid border check failed as D<0 " <<  std::endl;}
    else {
            mid_x1 = (-ell_c*sqr(ell_a*tan(theta_in)) + sqrt(mid_D)) / (sqr(ell_a*tan(theta_in))+ell_b*ell_b);
            mid_x2 = (-ell_c*sqr(ell_a*tan(theta_in)) - sqrt(mid_D)) / (sqr(ell_a*tan(theta_in))+ell_b*ell_b);
            std::cout << "mid check " << mid_x1 << " " << mid_x2 <<  std::endl;
            mid_x = std::max(mid_x1, mid_x2);
    }

    d_alpha = (ell_c+mid_x)/cos(theta_in);

    std::cout << "alpha " << d_alpha <<  std::endl; 
    ell_point_glob = p_focus + axtoF * d_alpha;
    std::cout << "ell_point_glob " << ell_point_glob[0] << " " << ell_point_glob[1] << " " << ell_point_glob[2] <<  std::endl; 


    ell_point_loc = ell_point_glob - ell_ax_cent;
    std::cout << "ell_point_loc " << ell_point_loc[0] << " " << ell_point_loc[1] << " " << ell_point_loc[2] <<  std::endl;

    new_ell_point_loc[0]=elax1*ell_point_loc;
    // new_ell_point_loc[1]=ell_point_loc[1];
    new_ell_point_loc[1]=0.0;
    
    new_ell_point_loc[2]=sin(acos(elax1*ell_point_loc/ell_point_loc.norm()))*ell_point_loc.norm();

    std::cout << "new_ell_point_loc " << new_ell_point_loc[0] << " " << new_ell_point_loc[1] << " " << new_ell_point_loc[2] <<  std::endl;
    
    axF11=-(axtoF % axPM2);
    axF11 /= axF11.norm();
    std::cout << "focal ax1 " << axF11[0] << " " << axF11[1] << " " << axF11[2] <<  std::endl;
    std::cout << "focal ax2 " << axPM2[0] << " " << axPM2[1] << " " << axPM2[2] <<  std::endl;

    firstfocus1.emplace_back(PlaneSurface(p_focus, axF11, -axPM2));
                             

    std::cout << "plotting 1st focus" << std::endl;

    // plot_field_on_given_surface(focal_zone, beam, focal_region, 100, "E", "focus_test", "step2");

    std::string filename_suffixST2 = "step2_strat_"+std::to_string(i);

    plot_field_on_given_surface(firstfocus1[i], firstrefs1[i], focal_region, 101, "E", "focus1", filename_suffixST2);

    std::cout << "focus plotted" << std::endl;

    std::cout <<"1st focus pass" << std::endl;

    //---------------------------------------------determine beam borders on elliptic mirror----------------------------------------
    
    low3d=parmirrors1[i].point(low2d);
    // std::cout << "par mirror low " << low3d[0] << " " << low3d[1] << " " << low3d[2] <<  std::endl;
    
   lowtoF=-(low3d-p_focus);
    lowtoF /= lowtoF.norm();
    // std::cout << "lowtoF " << lowtoF[0] << " " << lowtoF[1] << " " << lowtoF[2] <<  std::endl;

    lowtomid = acos(lowtoF*axtoF);
    // std::cout << "lowtomid " << lowtomid/M_PI*180 <<  std::endl;

    new_lowtomid = theta_in + lowtomid;
    // std::cout << "new_lowtomid " << new_lowtomid/M_PI*180 <<  std::endl;

    low_D = sqr(ell_b*ell_a)*(sqr(tan(new_lowtomid))*(ell_a*ell_a-ell_c*ell_c)+ell_b*ell_b);
    low_x = 0.0;
    
    if (low_D<0.0) {std::cout << "lower border check failed as D<0 " <<  std::endl;}
    else {
            low_x1 = (-ell_c*sqr(ell_a*tan(new_lowtomid)) + sqrt(low_D)) / (sqr(ell_a*tan(new_lowtomid))+ell_b*ell_b);
            low_x2 = (-ell_c*sqr(ell_a*tan(new_lowtomid)) - sqrt(low_D)) / (sqr(ell_a*tan(new_lowtomid))+ell_b*ell_b);
            low_x = std::max(low_x1, low_x2);
            if (low_x > new_ell_point_loc[0]) {std::cout << "lower border check failed as x doesn't exist " <<  std::endl;}
    }

    
    high3d=parmirrors1[i].point(high2d);
    // std::cout << "par mirror high " << high3d[0] << " " << high3d[1] << " " << high3d[2] <<  std::endl;
    
    hightoF=-(high3d-p_focus);
    hightoF /= hightoF.norm();
    // std::cout << "hightoF " << hightoF[0] << " " << hightoF[1] << " " << hightoF[2] <<  std::endl;

    hightomid = acos(hightoF*axtoF);
    // std::cout << "hightomid " << hightomid/M_PI*180 <<  std::endl;

    hightolow = acos(hightoF*lowtoF);
    std::cout << "hightolow " << hightolow/M_PI*180 <<  std::endl;

    new_hightomid = theta_in - hightomid;
    // std::cout << "new_hightomid " << new_hightomid/M_PI*180 <<  std::endl;


    high_D = sqr(ell_b*ell_a)*(sqr(tan(new_hightomid))*(ell_a*ell_a-ell_c*ell_c)+ell_b*ell_b);
    high_x = 0.0;
    
    if (high_D<0.0) {std::cout << "higher border check failed as D<0 " <<  std::endl;}
    else {
            high_x1 = (-ell_c*sqr(ell_a*tan(new_hightomid)) + sqrt(high_D)) / (sqr(ell_a*tan(new_hightomid))+ell_b*ell_b);
            high_x2 = (-ell_c*sqr(ell_a*tan(new_hightomid)) - sqrt(high_D)) / (sqr(ell_a*tan(new_hightomid))+ell_b*ell_b);
            high_x = std::max(high_x1, high_x2);
            if (high_x < new_ell_point_loc[0]) {std::cout << "higher border check failed as x doesn't exist " <<  std::endl;}
    }

    midtoF=MID3d-p_focus;
    Ftoellmid=p_focus-ell_point_glob;
    ell_mirr_length=beam_width*Ftoellmid.norm()/midtoF.norm();
    std::cout << "ell length " << ell_mirr_length << std::endl;
    //---------------------------------------------------------------------------------------------------------------------------------------------------

    ellregs1.emplace_back();
    // std::cout << "ellregs1.x_min " << ellregs1[i].x_min << std::endl;
    // std::cout << "ellregs1.y_max " << ellregs1[i].y_max << std::endl;
    
    // ell_beam_region.x_min = new_ell_point_loc[0] - d_alpha*0.2;
    // ell_beam_region.x_max = new_ell_point_loc[0] + d_alpha*0.15;
    // ell_beam_region.y_min = new_ell_point_loc[1] - d_alpha*0.15;
    // ell_beam_region.y_max = new_ell_point_loc[1] + d_alpha*0.15;

    std::cout << "highx " << high_x << std::endl;
    std::cout << "lowx " << low_x << std::endl;

    // Vector2D highlow2dcheck;
    // highlow2dcheck[0]=high_x;
    // highlow2dcheck[1]=-ell_mirr_length*0.5;
    // Position highlow3dcheck=ell_mirror.point(highlow2dcheck);
    // std::cout << "beam coords 1 " << highlow3dcheck[0] <<" "<< highlow3dcheck[1]  <<" "<< highlow3dcheck[2] << std::endl;
    // highlow2dcheck[0]=low_x;
    // highlow2dcheck[1]=-ell_mirr_length*0.5;
    // highlow3dcheck=ell_mirror.point(highlow2dcheck);
    // std::cout << "beam coords 2 " << highlow3dcheck[0] <<" "<< highlow3dcheck[1]  <<" "<< highlow3dcheck[2] << std::endl;
    // highlow2dcheck[0]=high_x;
    // highlow2dcheck[1]=ell_mirr_length*0.5;
    // highlow3dcheck=ell_mirror.point(highlow2dcheck);
    // std::cout << "beam coords 3 " << highlow3dcheck[0] <<" "<< highlow3dcheck[1]  <<" "<< highlow3dcheck[2] << std::endl;
    // highlow2dcheck[0]=low_x;
    // highlow2dcheck[1]=ell_mirr_length*0.5;
    // highlow3dcheck=ell_mirror.point(highlow2dcheck);
    // std::cout << "beam coords 4 " << highlow3dcheck[0] <<" "<< highlow3dcheck[1]  <<" "<< highlow3dcheck[2] << std::endl;

    ellregs1[i].x_min = low_x;
    ellregs1[i].x_max = high_x;
    ellregs1[i].y_min = -ell_mirr_length*0.5;
    ellregs1[i].y_max = ell_mirr_length*0.5;

    // std::cout << "ellregs1.x_min " << ellregs1[i].x_min << std::endl;
    // std::cout << "ellregs1.y_max " << ellregs1[i].y_max << std::endl;

    focus1_pass1.emplace_back(firstfocus1[i], firstrefs1[i], focal_region);

    // StrattonChuReflectionLowPrec focus1_pass_interp(focal_zone, Focus1Field, focal_region);



    // plot_field_on_given_surface(ell_mirror, beam, ell_beam_region, 100, "E", "mirror_ell_test", "step3");

    

    std::cout << "plotting ell mirror" << std::endl;

    std::string filename_suffixST3 = "step3_strat_"+std::to_string(i);

    plot_field_on_given_surface(ellmirrors1[i], focus1_pass1[i], ellregs1[i], 11, "E", "mirror_ell1_pass", filename_suffixST3);

    std::cout << "ell mirror plotted" << std::endl;


    std::cout <<"2nd reflect" << std::endl;


    secondrefs1.emplace_back(ellmirrors1[i], focus1_pass1[i], ellregs1[i]);

    std::cout << "second_focus " << second_focus[0] << " " << second_focus[1] << " " << second_focus[2] <<  std::endl;


    // Vector axtoF2=second_focus-ell_point_glob;
    // axtoF2 /= axtoF2.norm();
    // std::cout << "axtoF2 " << axtoF2[0] << " " << axtoF2[1] << " " << axtoF2[2] <<  std::endl;

    // Vector axF21=-(axtoF2 % axPM2);
    // axF21 /= axF21.norm();
    // std::cout << "focal ax1 " << axF21[0] << " " << axF21[1] << " " << axF21[2] <<  std::endl;

    // PlaneSurface focal_zone2(second_focus, axF21, axPM2);


    // PlaneSurface second_focal_plane(second_focus, elax1, elax3);
    // PlaneSurface first_focal_plane(p_focus, elax1, elax3);


    // plot_field_on_given_surface(second_focal_plane, beam, focal_region, 100, "E", "focus2_test", "step4");
    // plot_field_on_given_surface(first_focal_plane, beam, focal_region, 100, "E", "focus1_test", "step4");
    // plot_field_on_given_surface(focal_zone2, beam, focal_region, 100, "E", "focus2_norm_test", "step4");

    // StrattonChuReflection ellfocus1_pass(first_focal_plane, reflect1, focal_region);

    // plot_field_on_given_surface(first_focal_plane, reflect1, focal_region, 100, "E", "NEWfocus1_test", "step4");

    std::cout << "plotting ell mirror second time" << std::endl;

    // plot_field_on_given_surface(ell_mirror, ellfocus1_pass, ell_beam_region, 10, "E", "NEWmirror_ell_pass", "step3_strat");
    // StrattonChuReflection ellreflect2(ell_mirror, ellfocus1_pass, ell_beam_region);

    std::cout << "plotting 2nd focus" << std::endl;
    // plot_field_on_given_surface(second_focal_plane, ellreflect2, focal_region, second_focus_points, "E", "second_focus", "step4");

    // plot_field_on_given_surface(second_focal_plane, reflect2, focal_region, second_focus_points, "E", "second_focus", "step4");

        

    // plot_field_on_given_surface_ell(focal_zone2, reflect2, focal_regionELL2, second_focus_points, "E", "second_focus", "step4");
    // plot_field_on_given_surface(second_focal_plane, reflectF, focal_region, 111, "E", "second_focus", "stepTT");

    std::cout <<"2nd focus plotted" << std::endl;
    // if (i==0) {
    //     Vector axtoF2=second_focus-ell_point_glob;
    //     axtoF2 /= axtoF2.norm();
    //     std::cout << "axtoF2 " << axtoF2[0] << " " << axtoF2[1] << " " << axtoF2[2] <<  std::endl;

    //     Vector axF21=-(axtoF2 % axPM2);
    //     axF21 /= axF21.norm();
    //     std::cout << "focal ax1 " << axF21[0] << " " << axF21[1] << " " << axF21[2] <<  std::endl;

    //     PlaneSurface focal_zone2(second_focus, axF21, axPM2);
    //     std::cout << "plotting 2nd focus" << std::endl;
        
    //     plot_field_on_given_surface(focal_zone2, incbeams1[0], focal_regionELL2, 100, "E", "focus2_norm_test", "step4");
    //     plot_field_on_given_surface_ell(focal_zone2, secondrefs1[0], focal_regionELL2, 11, "E", "second_focus", "step4");
    //     // plot_field_on_given_surface(second_focal_plane, reflectF, focal_region, 111, "E", "second_focus", "stepTT");

    //     std::cout <<"2nd focus plotted" << std::endl;
    // }


    }



    Vector2D ell_high2d;
    ell_high2d[0]=high_x;
    ell_high2d[1]=0;
    Position ell_high3d=ellmirrors1[5].point(ell_high2d);

    Vector2D ell_low2d;
    ell_low2d[0]=low_x;
    ell_low2d[1]=0;
    Position ell_low3d=ellmirrors1[5].point(ell_low2d);

    Vector ell_hightof2=second_focus-ell_high3d;
    Vector ell_lowtof2=second_focus-ell_low3d;
    Vector ell_mirr_height=ell_high3d-ell_low3d;

    double emh=ell_mirr_height.norm();
    double emhtf=ell_hightof2.norm();
    double emltf=ell_lowtof2.norm();


    std::cout << "ell_mirr high to F2 " << emhtf <<  std::endl;
    std::cout << "ell_mirr low to F2 " << emltf <<  std::endl;
    std::cout << "ell_mirr height " << emh <<  std::endl;

    double shozhd_angle=acos((emhtf*emhtf+emltf*emltf-emh*emh)/(2*emhtf*emltf));
    std::cout << "ugol shozhdeniya " << shozhd_angle/M_PI*180 <<  std::endl;


    ell_hightof2 /= ell_hightof2.norm();
    ell_lowtof2 /= ell_lowtof2.norm();
    double clever_shozhd_angle=acos(ell_hightof2*ell_lowtof2);
    std::cout << "ugol shozhdeniya umom " << clever_shozhd_angle/M_PI*180 <<  std::endl;

    axFm1 /= -axFm1.norm();
    double psi_off=acos(axtoF*axFm1);
    double real_focal=2*F/(1+axtoF*axFm1);
    std::cout << "psi_off " << psi_off/M_PI*180 <<  std::endl;
    std::cout << "alpha " << alpha/M_PI*180 <<  std::endl;
    std::cout << "real_focal " << real_focal <<  std::endl;



    for (size_t i=0;i<beams_number;i++){
    // size_t i=0;

    axFm2={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
    ax1m2={-1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), -1.0*cos(gamma_med)};
    ax2m2=(ax1m2 % axFm2);
    ax2m2 /= ax2m2.norm();

    parmirrors2.emplace_back(ParabolicSurface(
            {-(F*cos(gamma_med)+addlength)*cos(2*M_PI/beams_number*i+add_angle), -(F*cos(gamma_med)+addlength)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)-addheight},
            ax1m2, ax2m2, 4.0*F, 4.0*F
        ));
    ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), 1.0*cos(gamma_med)};

    incbeams2.emplace_back(ParallelBeamAlpha(lambda, {-addlength*cos(2*M_PI/beams_number*i+add_angle), -addlength*sin(2*M_PI/beams_number*i+add_angle), -addheight},/*new_r02,*/ ax1m1, ax2m2,
                                //  [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 + h) + sqr(x2)) / (2 * sqr(sigma)) ); },
                                // [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
                                // [h, gauss_A, sigma, zern=ZernikeAberrations<Order>(sigma_poly[sig_ind], beam_radius*sqrt(2), i+beams_number)](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*zern(x1+h, x2)); },
                                // [h, gauss_A, sigma, i, beams_number, &beam_aberr, zern_ind_n, zern_ind_m](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i+beams_number]*ZernikeCartesian(zern_ind_n,zern_ind_m,x1+h,x2, beam_radius*sqrt(2))); },
                                [h, gauss_A, sigma, i, beams_number, &beam_aberr, zern_ind_j](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ) * exp(std::complex<double>(0.0, 1.0)*beam_aberr[i+beams_number]*ZernikeCartesianOnSquare1(zern_ind_j,x1+h,x2, beam_radius));},
                                //[h](double x1, double x2) { return sqr(x1 + h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
                                // [h](double x1, double x2) { return smoothed(sqrt(sqr(x1 + h) + sqr(x2)), beam_radius); },
                                [](double, double) { return 0.0; }));
    

    p_focus = {-addlength*cos(2*M_PI/beams_number*i+add_angle), -addlength*sin(2*M_PI/beams_number*i+add_angle), -addheight};
    Vector elax3 = -ax2m2;
    Vector axPM2 = ax2m2;


    elax1=main_focus-p_focus;
    elax1_copy=main_focus-p_focus;
    elax1 /= elax1.norm();


    std::cout << "elax1 " << elax1[0] << " " << elax1[1] << " " << elax1[2] <<  std::endl;
    std::cout << "elax3 " << elax3[0] << " " << elax3[1] << " " << elax3[2] <<  std::endl;

    ell_c=elax1_copy.norm()/2;
    ell_ax_cent = p_focus + elax1 * ell_c;
    std::cout << "ell_ax_cent " << ell_ax_cent[0] << " " << ell_ax_cent[1] << " " << ell_ax_cent[2] <<  std::endl;
    

    ellmirrors2.emplace_back(EllipticSurface(ell_ax_cent, elax1, elax3, ell_a, ell_b, ell_b)); 


  
    std::cout << "plotting parab mirror" << std::endl;

    std::string filename_suffixST1 = "step1_"+std::to_string(i+beams_number);

    plot_field_on_given_surface(parmirrors2[i], incbeams2[i], region1g, 111, "E", "mirror_par2", filename_suffixST1);

    std::cout << "parab mirror plotted" << std::endl;
    
    firstrefs2.emplace_back(parmirrors2[i], incbeams2[i], region1g);

    
    MID3d=parmirrors2[i].point(MID2d);
    std::cout << "mirror mid " << MID3d[0] << " " << MID3d[1] << " " << MID3d[2] <<  std::endl;

    std::cout << "p_focus " << p_focus[0] << " " << p_focus[1] << " " << p_focus[2] <<  std::endl;

    
    axtoF=-(MID3d-p_focus);
    axtoF /= axtoF.norm();
    std::cout << "axtoF " << axtoF[0] << " " << axtoF[1] << " " << axtoF[2] <<  std::endl;

    theta_in = acos(elax1*axtoF);

    std::cout << "theta_in "<< theta_in/M_PI*180 <<  std::endl;

    mid_D = sqr(ell_b*ell_a)*(sqr(tan(theta_in))*(ell_a*ell_a-ell_c*ell_c)+ell_b*ell_b);
    mid_x = 0.0;
    
    if (mid_D<0.0) {std::cout << "mid border check failed as D<0 " <<  std::endl;}
    else {
            mid_x1 = (-ell_c*sqr(ell_a*tan(theta_in)) + sqrt(mid_D)) / (sqr(ell_a*tan(theta_in))+ell_b*ell_b);
            mid_x2 = (-ell_c*sqr(ell_a*tan(theta_in)) - sqrt(mid_D)) / (sqr(ell_a*tan(theta_in))+ell_b*ell_b);
            std::cout << "mid check " << mid_x1 << " " << mid_x2 <<  std::endl;
            mid_x = std::max(mid_x1, mid_x2);
    }

    d_alpha = (ell_c+mid_x)/cos(theta_in);

    std::cout << "alpha " << d_alpha <<  std::endl; 
    ell_point_glob = p_focus + axtoF * d_alpha;
    std::cout << "ell_point_glob " << ell_point_glob[0] << " " << ell_point_glob[1] << " " << ell_point_glob[2] <<  std::endl; 


    ell_point_loc = ell_point_glob - ell_ax_cent;
    std::cout << "ell_point_loc " << ell_point_loc[0] << " " << ell_point_loc[1] << " " << ell_point_loc[2] <<  std::endl;

    new_ell_point_loc[0]=elax1*ell_point_loc;
    // new_ell_point_loc[1]=ell_point_loc[1];
    new_ell_point_loc[1]=0.0;
    
    new_ell_point_loc[2]=sin(acos(elax1*ell_point_loc/ell_point_loc.norm()))*ell_point_loc.norm();

    std::cout << "new_ell_point_loc " << new_ell_point_loc[0] << " " << new_ell_point_loc[1] << " " << new_ell_point_loc[2] <<  std::endl;
    
    axF11=-(axtoF % axPM2);
    axF11 /= axF11.norm();
    std::cout << "focal ax1 " << axF11[0] << " " << axF11[1] << " " << axF11[2] <<  std::endl;
    std::cout << "focal ax2 " << axPM2[0] << " " << axPM2[1] << " " << axPM2[2] <<  std::endl;

    firstfocus2.emplace_back(PlaneSurface(p_focus, axF11, -axPM2));
                             

    std::cout << "plotting 1st focus" << std::endl;

    // plot_field_on_given_surface(focal_zone, beam, focal_region, 100, "E", "focus_test", "step2");

    std::string filename_suffixST2 = "step2_strat_"+std::to_string(i+beams_number);

    plot_field_on_given_surface(firstfocus2[i], firstrefs2[i], focal_region, 101, "E", "focus2", filename_suffixST2);

    std::cout << "focus plotted" << std::endl;

    std::cout <<"1st focus pass" << std::endl;

    //---------------------------------------------determine beam borders on elliptic mirror----------------------------------------
    
    low3d=parmirrors2[i].point(low2d);
    // std::cout << "par mirror low " << low3d[0] << " " << low3d[1] << " " << low3d[2] <<  std::endl;
    
   lowtoF=-(low3d-p_focus);
    lowtoF /= lowtoF.norm();
    // std::cout << "lowtoF " << lowtoF[0] << " " << lowtoF[1] << " " << lowtoF[2] <<  std::endl;

    lowtomid = acos(lowtoF*axtoF);
    // std::cout << "lowtomid " << lowtomid/M_PI*180 <<  std::endl;

    new_lowtomid = theta_in + lowtomid;
    // std::cout << "new_lowtomid " << new_lowtomid/M_PI*180 <<  std::endl;

    low_D = sqr(ell_b*ell_a)*(sqr(tan(new_lowtomid))*(ell_a*ell_a-ell_c*ell_c)+ell_b*ell_b);
    low_x = 0.0;
    
    if (low_D<0.0) {std::cout << "lower border check failed as D<0 " <<  std::endl;}
    else {
            low_x1 = (-ell_c*sqr(ell_a*tan(new_lowtomid)) + sqrt(low_D)) / (sqr(ell_a*tan(new_lowtomid))+ell_b*ell_b);
            low_x2 = (-ell_c*sqr(ell_a*tan(new_lowtomid)) - sqrt(low_D)) / (sqr(ell_a*tan(new_lowtomid))+ell_b*ell_b);
            low_x = std::max(low_x1, low_x2);
            if (low_x > new_ell_point_loc[0]) {std::cout << "lower border check failed as x doesn't exist " <<  std::endl;}
    }

    
    high3d=parmirrors2[i].point(high2d);
    // std::cout << "par mirror high " << high3d[0] << " " << high3d[1] << " " << high3d[2] <<  std::endl;
    
    hightoF=-(high3d-p_focus);
    hightoF /= hightoF.norm();
    // std::cout << "hightoF " << hightoF[0] << " " << hightoF[1] << " " << hightoF[2] <<  std::endl;

    hightomid = acos(hightoF*axtoF);
    // std::cout << "hightomid " << hightomid/M_PI*180 <<  std::endl;

    hightolow = acos(hightoF*lowtoF);
    std::cout << "hightolow " << hightolow/M_PI*180 <<  std::endl;

    new_hightomid = theta_in - hightomid;
    // std::cout << "new_hightomid " << new_hightomid/M_PI*180 <<  std::endl;


    high_D = sqr(ell_b*ell_a)*(sqr(tan(new_hightomid))*(ell_a*ell_a-ell_c*ell_c)+ell_b*ell_b);
    high_x = 0.0;
    
    if (high_D<0.0) {std::cout << "higher border check failed as D<0 " <<  std::endl;}
    else {
            high_x1 = (-ell_c*sqr(ell_a*tan(new_hightomid)) + sqrt(high_D)) / (sqr(ell_a*tan(new_hightomid))+ell_b*ell_b);
            high_x2 = (-ell_c*sqr(ell_a*tan(new_hightomid)) - sqrt(high_D)) / (sqr(ell_a*tan(new_hightomid))+ell_b*ell_b);
            high_x = std::max(high_x1, high_x2);
            if (high_x < new_ell_point_loc[0]) {std::cout << "higher border check failed as x doesn't exist " <<  std::endl;}
    }

    midtoF=MID3d-p_focus;
    Ftoellmid=p_focus-ell_point_glob;
    ell_mirr_length=beam_width*Ftoellmid.norm()/midtoF.norm();
    std::cout << "ell length " << ell_mirr_length << std::endl;
    //---------------------------------------------------------------------------------------------------------------------------------------------------

    ellregs2.emplace_back();

    // ell_beam_region.x_min = new_ell_point_loc[0] - d_alpha*0.2;
    // ell_beam_region.x_max = new_ell_point_loc[0] + d_alpha*0.15;
    // ell_beam_region.y_min = new_ell_point_loc[1] - d_alpha*0.15;
    // ell_beam_region.y_max = new_ell_point_loc[1] + d_alpha*0.15;

    std::cout << "highx " << high_x << std::endl;
    std::cout << "lowx " << low_x << std::endl;

    // Vector2D highlow2dcheck;
    // highlow2dcheck[0]=high_x;
    // highlow2dcheck[1]=-ell_mirr_length*0.5;
    // Position highlow3dcheck=ell_mirror.point(highlow2dcheck);
    // std::cout << "beam coords 1 " << highlow3dcheck[0] <<" "<< highlow3dcheck[1]  <<" "<< highlow3dcheck[2] << std::endl;
    // highlow2dcheck[0]=low_x;
    // highlow2dcheck[1]=-ell_mirr_length*0.5;
    // highlow3dcheck=ell_mirror.point(highlow2dcheck);
    // std::cout << "beam coords 2 " << highlow3dcheck[0] <<" "<< highlow3dcheck[1]  <<" "<< highlow3dcheck[2] << std::endl;
    // highlow2dcheck[0]=high_x;
    // highlow2dcheck[1]=ell_mirr_length*0.5;
    // highlow3dcheck=ell_mirror.point(highlow2dcheck);
    // std::cout << "beam coords 3 " << highlow3dcheck[0] <<" "<< highlow3dcheck[1]  <<" "<< highlow3dcheck[2] << std::endl;
    // highlow2dcheck[0]=low_x;
    // highlow2dcheck[1]=ell_mirr_length*0.5;
    // highlow3dcheck=ell_mirror.point(highlow2dcheck);
    // std::cout << "beam coords 4 " << highlow3dcheck[0] <<" "<< highlow3dcheck[1]  <<" "<< highlow3dcheck[2] << std::endl;

    ellregs2[i].x_min = low_x;
    ellregs2[i].x_max = high_x;
    ellregs2[i].y_min = -ell_mirr_length*0.5;
    ellregs2[i].y_max = ell_mirr_length*0.5;

    focus1_pass2.emplace_back(firstfocus2[i], firstrefs2[i], focal_region);

    // StrattonChuReflectionLowPrec focus1_pass_interp(focal_zone, Focus1Field, focal_region);



    // plot_field_on_given_surface(ell_mirror, beam, ell_beam_region, 100, "E", "mirror_ell_test", "step3");

    

    std::cout << "plotting ell mirror" << std::endl;

    std::string filename_suffixST3 = "step3_strat_"+std::to_string(i+beams_number);

    plot_field_on_given_surface(ellmirrors2[i], focus1_pass2[i], ellregs2[i], 11, "E", "mirror_ell2_pass", filename_suffixST3);

    std::cout << "ell mirror plotted" << std::endl;

    std::cout <<"2nd reflect" << std::endl;


    secondrefs2.emplace_back(ellmirrors2[i], focus1_pass2[i], ellregs2[i]);

    std::cout << "second_focus " << second_focus[0] << " " << second_focus[1] << " " << second_focus[2] <<  std::endl;


    // Vector axtoF2=second_focus-ell_point_glob;
    // axtoF2 /= axtoF2.norm();
    // std::cout << "axtoF2 " << axtoF2[0] << " " << axtoF2[1] << " " << axtoF2[2] <<  std::endl;

    // Vector axF21=-(axtoF2 % axPM2);
    // axF21 /= axF21.norm();
    // std::cout << "focal ax1 " << axF21[0] << " " << axF21[1] << " " << axF21[2] <<  std::endl;

    // PlaneSurface focal_zone2(second_focus, axF21, axPM2);


    // PlaneSurface second_focal_plane(second_focus, elax1, elax3);
    // PlaneSurface first_focal_plane(p_focus, elax1, elax3);


    // plot_field_on_given_surface(second_focal_plane, beam, focal_region, 100, "E", "focus2_test", "step4");
    // plot_field_on_given_surface(first_focal_plane, beam, focal_region, 100, "E", "focus1_test", "step4");
    // plot_field_on_given_surface(focal_zone2, beam, focal_region, 100, "E", "focus2_norm_test", "step4");

    // StrattonChuReflection ellfocus1_pass(first_focal_plane, reflect1, focal_region);

    // plot_field_on_given_surface(first_focal_plane, reflect1, focal_region, 100, "E", "NEWfocus1_test", "step4");

    std::cout << "plotting ell mirror second time" << std::endl;

    // plot_field_on_given_surface(ell_mirror, ellfocus1_pass, ell_beam_region, 10, "E", "NEWmirror_ell_pass", "step3_strat");
    // StrattonChuReflection ellreflect2(ell_mirror, ellfocus1_pass, ell_beam_region);

    std::cout << "plotting 2nd focus" << std::endl;
    // plot_field_on_given_surface(second_focal_plane, ellreflect2, focal_region, second_focus_points, "E", "second_focus", "step4");

    // plot_field_on_given_surface(second_focal_plane, reflect2, focal_region, second_focus_points, "E", "second_focus", "step4");

        

    // plot_field_on_given_surface_ell(focal_zone2, reflect2, focal_regionELL2, second_focus_points, "E", "second_focus", "step4");
    // plot_field_on_given_surface(second_focal_plane, reflectF, focal_region, 111, "E", "second_focus", "stepTT");

    std::cout <<"2nd focus plotted" << std::endl;
    // std::cout << i << std::endl;
    }

    



//     // // // ----------------------------------------------------------------Building reflections-------------------------------------------------------------------
    

    // std::vector<StrattonChuReflectionLowPrec> refs;
    // refs.reserve(1);
    // refs.push_back(secondrefs1[0]);
    // refs.emplace_back(mirror, beam, region1g);

    std::vector<StrattonChuReflectionLowPrec> refs;
    //double belt
    refs.reserve(2*beams_number);
    for (size_t i=0;i<2*beams_number;i++){
        if (i<beams_number){
            refs.push_back(secondrefs1[i]);
        }
        else {
            refs.push_back(secondrefs2[i-beams_number]);
        }        
    }

//     // //single belt
//     // refs.reserve(beams_number);
//     // for (size_t i=0;i<beams_number;i++){
//     //     refs.emplace_back(mirrors1[i], beams1[i], region1);              
//     // }

    // std::vector<StrattonChuReflectionLowPrec> refs;
    // refs.reserve(1);
    // refs.emplace_back(ell_mirror, focus1_pass, ell_beam_region);
    // std::string nameEL="second_focus_21_points";
    // FieldsInRegion <21> fields_long(focal_zone2, focal_regionELL2, nameEL, refs, 1);   
 

    // PrlppdVolume focal_volume( p_focus, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
    // PrlppdVolume focal_volume( main_focus, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
    // VolumeRegion focal_area;
    // double xcent=0.0;
    // // double xcent=8 * lambda;
    // // double xcent=-8 * lambda;
    // double ycent=0.0;
    // // double ycent=8 * lambda;
    // // double ycent=-8 * lambda;
    // double zcent=0.0;
    // // double zcent=8 * lambda;
    // // double zcent=-8 * lambda;
    // focal_area.x_min = xcent - 2 * lambda;
    // focal_area.x_max = xcent + 2 * lambda;
    // focal_area.y_min = ycent - 2 * lambda;
    // focal_area.y_max = ycent + 2 * lambda;
    // focal_area.z_min = zcent - 2 * lambda;
    // focal_area.z_max = zcent + 2 * lambda;


    PlaneSurface focal_plane_long( main_focus, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
    PlaneSurface focal_plane_trans( main_focus, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0});
    SurfaceRegion focal_reg;
    double xcent=0.0;
    // double xcent=8 * lambda;
    // double xcent=-8 * lambda;
    double ycent=0.0;
    // double ycent=8 * lambda;
    // double ycent=-8 * lambda;
    focal_reg.x_min = xcent - 1.5 * lambda;
    focal_reg.x_max = xcent + 1.5 * lambda;
    focal_reg.y_min = ycent - 1.5 * lambda;
    focal_reg.y_max = ycent + 1.5 * lambda;

//     // std::cout << "meow" << std::endl;
//     // //double belt

//------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // std::string nameE=std::string(GET_VARIABLE_NAME(focal_area))+"_"+std::to_string(2*beams_number)+"-beams_new_config";
    // // std::string nameE=std::string(GET_VARIABLE_NAME(focal_area))+"_"+std::to_string(1)+"-beam_in_12config";
    // // // std::cout << nameE << std::endl;
    // FieldsIn3dRegion <3> fields(focal_volume, focal_area, nameE, refs, 2*beams_number);
    // // FieldsIn3dRegion <focal_points> fields(focal_volume, focal_area, nameE, refs, 2*beams_number);
//     // FieldsIn3dRegion <focal_points> fields(focal_volume, focal_area, nameE, refs, 1);
    // FieldsIn3dRegion <5> fields(focal_volume, focal_area, nameE, refs, 1);
//     std::cout <<"all area ready" << std::endl;

    // std::string name=std::string(GET_VARIABLE_NAME(focal_plane_long))+"_"+std::to_string(2*beams_number)+"-beams_new_config";
    std::string name=std::string(GET_VARIABLE_NAME(focal_plane_long))+"_"+std::to_string(2*beams_number)+"-beams_new_config_with_aberrs";
    // std::string name=std::string(GET_VARIABLE_NAME(focal_zone))+"_"+std::to_string(1)+"-beam_in_12config";
    // std::string name="first_focus";
    FieldsInRegion <7> fields_long(focal_plane_long, focal_reg, name, refs, 2*beams_number);
    // FieldsInRegion <14> fields_long(focal_zone, focal_region, name, refs, 1);
    // FieldsInRegion <focal_points> fields_long(focal_plane, focal_region, name, refs, 1);          
    std::cout <<"fp ready" << std::endl;
    // name=std::string(GET_VARIABLE_NAME(focal_plane_trans))+"_"+std::to_string(2*beams_number)+"-beams_new_config";
    name=std::string(GET_VARIABLE_NAME(focal_plane_trans))+"_"+std::to_string(2*beams_number)+"-beams_new_config_with_aberrs";
    FieldsInRegion <7> fields_trans(focal_plane_trans, focal_reg, name, refs, 2*beams_number);
    std::cout <<"fp trans ready" << std::endl;

// //     // // //single belt
// //     // // std::string name=std::string(GET_VARIABLE_NAME(focal_region))+"_"+std::to_string(beams_number)+"-beams_single";
// //     // // FieldsInRegion <focal_points> fields_long(focal_plane, focal_region, name, refs, beams_number);    
// //     // // std::cout <<"fp ready" << std::endl;
// //     // // name=std::string(GET_VARIABLE_NAME(focal_region_transversal))+"_"+std::to_string(beams_number)+"-beams_single";
// //     // // FieldsInRegion <focal_points_transversal> fields_trans(focal_plane_transversal, focal_region_transversal, name, refs, beams_number);
// //     // // std::cout <<"fp trans ready" << std::endl;

// //     // // std::vector<double> phase_lims(6);
// //     // // phase_lims[0]=M_PI/3;
// //     // // phase_lims[1]=M_PI/4;
// //     // // phase_lims[2]=M_PI/6;
// //     // // phase_lims[3]=0;
// //     // // phase_lims[4]=M_PI/2;
// //     // // phase_lims[5]=M_PI;

 


    std::vector<double> phase_lims(5);
    phase_lims[0]=0;
    phase_lims[1]=M_PI/3;
    phase_lims[2]=M_PI/2;
    phase_lims[3]=M_PI;
    phase_lims[4]=2*M_PI;

    std::vector<double> phases;
    phases.reserve(2*beams_number);
//     // phases.reserve(beams_number);

    std::random_device dev;
    std::mt19937 rng(dev());

    double testIntens[5][10];
    double trueint;

    std::string intens_file = "intens_for_"+std::to_string(2*beams_number)+"_beams.txt";
    // std::string intens_file = "intens_for_1beam_elliptic.txt";
    // std::string intens_file = "intens_for_"+std::to_string(beams_number)+"_beams_single.txt";
    std::ofstream outr(intens_file);


//     // std::string filename_suffix1;
//     // for (size_t i=0;i<2*beams_number;i++){
//     //     if (i<beams_number){
//     //         filename_suffix1 = "_beam"+std::to_string(i+1)+"_side1";
//     //     }
//     //     else {
//     //         filename_suffix1 = "_beam"+std::to_string(i+1-beams_number)+"_side2";
//     //     }    
//     //     // plot_field_on_given_surface(planes[i], refs[i], focal_region, focal_points, "E1", "longitudinal_E1", filename_suffix1);
//     //     plot_field_on_given_surface_3darr(planes[i], focal_region, fields_long, 0, focal_points, "E1", "longitudinal_E1", filename_suffix1);
//     //     // plot_field_on_given_surface_3darr(focal_plane, focal_region, fields_long, i, focal_points, "E1", "longitudinal_E1", filename_suffix1);
//     //     // plot_field_on_given_surface_3darr(focal_plane_transversal, focal_region_transversal, fields_trans, i, focal_points_transversal, "E1", "transversal_E1", filename_suffix1);
//     //     std::cout << "meow" << std::endl;
//     // }
    outr << " phi = " << phi/M_PI*180 << '\n';
    

    for (size_t phase_iter=0;phase_iter<1;phase_iter++){
        std::uniform_real_distribution<double> dist(-phase_lims[phase_iter]/2,phase_lims[phase_iter]/2);
        if (phase_iter==4 || phase_iter==0){
            std::uniform_real_distribution<double> dist(-phase_lims[phase_iter]/2,phase_lims[phase_iter]/2);
            std::cout <<"max = " <<  phase_lims[phase_iter] << std::endl;
            outr << "for max = " <<  phase_lims[phase_iter] << '\n';
        } else {
            std::normal_distribution<double> dist{0, phase_lims[phase_iter]/2};
            std::cout <<"sigma = " <<  phase_lims[phase_iter]/2 << std::endl;
            outr << "for sigma = " <<  phase_lims[phase_iter]/2 << '\n';
        }
        
        size_t test_iter=0;
        // for (size_t test_iter=0; test_iter<10;test_iter++){
            for (size_t i=0;i<2*beams_number;i++){
            // for (size_t i=0;i<beams_number;i++){
                // std::cout <<"file" << i  << std::endl;
                phases[i]=dist(rng);
                // phases[i]=0;
                // std::cout << phases[i] << std::endl;
            }
            std::string filename_suffix = "_"+std::to_string(2*beams_number)+"beams_test"+std::to_string(phase_iter)+"_"+std::to_string(test_iter);
            // std::string filename_suffix = "_"+std::to_string(beams_number)+"beams_single_test"+std::to_string(phase_iter)+"_"+std::to_string(test_iter);

            // fields_long.constructSummaryMono(phases);
            // plot_field_on_given_surface_3darr(focal_plane, focal_region, fields_long, focal_points, "E1", "longitudinal_E1", filename_suffix);
            // fields_long.constructSummaryMono(phases);
            // plot_field_on_given_surface_3darr(focal_zone, focal_region, fields_long, 14, "E1", "longitudinal_E1", filename_suffix);

            fields_long.constructSummaryMono(phases);
            plot_field_on_given_surface_3darr(focal_plane_long, focal_reg, fields_long, 7, "E1", "longitudinal_E1", filename_suffix);
            fields_trans.constructSummaryMono(phases);
            plot_field_on_given_surface_3darr(focal_plane_trans, focal_reg, fields_trans, 7, "E1", "transversal_E1", filename_suffix);
            


            // plot_field_on_given_surface_3darr(focal_zone2, focal_regionELL2, fields_long, 21, "E_sf", "longitudinal_E1", filename_suffix);
            // plot_field_on_given_surface_3darr(focal_plane, focal_region, fields_long, focal_points, "E1", "longitudinal_E1", filename_suffix);
            // fields_trans.constructSummaryMono(phases);
            // plot_field_on_given_surface_3darr(focal_plane_transversal, focal_region_transversal, fields_trans, focal_points_transversal, "E1", "transversal_E1", filename_suffix);
            // fields.constructSummaryMono(phases);
            // plot_field_in_given_area_3darr(focal_volume, focal_area, fields, 3, "E1", "field", filename_suffix);

            fields_long.findSummaryMax();
            fields_long.findSummaryMaxExact();
            if (fields_long.checkMax()){
                std::cout <<"ok"<< std::endl;
            }
            else{
                std::cout <<"not ok"<< std::endl;
                trueint=fields_long.calculateTrueMaxIntensitySummary();
                std::cout <<"test "<< test_iter <<" max intens exact = " << trueint << std::endl;
            }
            testIntens[phase_iter][test_iter]=fields_long.calculateMaxIntensitySummary();
            std::cout <<"test "<< test_iter <<" max intens = " << testIntens[phase_iter][test_iter] << std::endl;
            Vector2D MaxPos=fields_long.MaxIntensityPoint();            
            outr << testIntens[phase_iter][test_iter] << "   in point [" << MaxPos[0] << ", " << MaxPos[1] <<"]" <<'\n';

            // fields.findSummaryMax();
            // fields.findSummaryMaxExact();
            // if (fields.checkMax()){
            //     std::cout <<"ok"<< std::endl;
            // }
            // else{
            //     std::cout <<"not ok"<< std::endl;
            //     trueint=fields.calculateTrueMaxIntensitySummary();
            //     std::cout <<"test "<< test_iter <<" max intens exact = " << trueint << std::endl;
            // }
            // testIntens[phase_iter][test_iter]=fields.calculateMaxIntensitySummary();
            // std::cout <<"test "<< test_iter <<" max intens = " << testIntens[phase_iter][test_iter] << std::endl;
            // Vector MaxPos=fields.MaxIntensityPoint();
            
            // outr << testIntens[phase_iter][test_iter] << "   in point [" << MaxPos[0] << ", " << MaxPos[1] << ", " << MaxPos[2] <<"]" <<'\n';

        // }
    }
    
//     outr.close();
  
    return 0;

}
