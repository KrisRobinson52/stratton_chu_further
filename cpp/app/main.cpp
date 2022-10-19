#include "stratton-chu/plane-surface.hpp"
#include "stratton-chu/parabolic-surface.hpp"
#include "stratton-chu/distorted-surface.hpp"
#include "stratton-chu/parallel-beam.hpp"
#include "stratton-chu/stratton-chu-field.hpp"
#include "stratton-chu/utils.hpp"
#include "stratton-chu/spec-inverse-fourier.hpp"

#include "stratton-chu/csv-saver.hpp"
#include "stratton-chu/vtk-saver.hpp"

#include <iostream>
#include "thread_pool.hpp"

#include <cmath>
#include <fftw3.h>

#include <fstream>

#include <vector>

#include "constants.hpp"

const double lambda = 0.000091; // cm = 910 nm
const double beam_radius = 10; // cm
const double beam_width = 2 * beam_radius; // cm



std::vector<SurfaceDistortionLegendre::DistortionPolinom> make_distortion_polynoms(double amp)
{
    std::vector<double> angles = {
        5.285787050983075,
        3.7315163082383522,
        0.4840372224686504,
        5.146130391061815,
        3.127904939599717,
        0.3754736274237447,
        5.7929281053760855,
        1.8623574557673694,
        3.902414951488807,
        4.844839876845448
    };

    std::vector<SurfaceDistortionLegendre::DistortionPolinom> result;

    for (size_t i = 0; i < angles.size(); i++)
    {
        result.push_back(SurfaceDistortionLegendre::DistortionPolinom(amp, angles[i], i*3+1));
    }

    return result;
}

void plot_field_on_given_surface(
        const ISurface& surface,
        const IField& field,
        const SurfaceRegion& region,
        int n_points,
        const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
{
    VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
    VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());

    // В тау-тау-эн;
    VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

    function task{[&region, &field, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m, n_points, i = 0, j = 0](){}};
    using TaskQueue = ConcurentQueue<decltype(task)>;
    using ThreadPool = ThreadPoolImpl<TaskQueue>;

    TaskQueue tasks;
    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            tasks.push([&region, &field, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m, n_points, i, j](){
                Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                Position p = surface.point(xy);
                FieldValue field_value = field.get(p);
                //std::cout << "point ready" << std::endl;

                vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
                vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
                vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

                FieldValue field_rotated;
                field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
                field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

                vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
                vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
            });
        }
    }
    ThreadPool threads{tasks};
    threads.run(worker<TaskQueue::value_type>);

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

    vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
}

void plot_field_on_given_surface_with_time(
        const ISurface& surface,
        const SpecIFT& field,
        const SurfaceRegion& region,
        int n_points,
        double timing,
        const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
{
    VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
    VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());

    // В тау-тау-эн;
    VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

    function task{[&region, &field, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m, n_points, timing, i = 0, j = 0](){}};
    using TaskQueue = ConcurentQueue<decltype(task)>;
    using ThreadPool = ThreadPoolImpl<TaskQueue>;

    TaskQueue tasks;
    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            tasks.push([&region, &field, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m, n_points, timing, i, j](){
                Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                Position p = surface.point(xy);
                FieldValue field_value = field.get(p, timing);
                //std::cout << "point ready "<< i << "," << j << std::endl;

                vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
                vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
                vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

                FieldValue field_rotated;
                field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
                field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

                vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
                vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
            });
        }
    }
    ThreadPool threads{tasks};
    threads.run(worker<TaskQueue::value_type>);

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

    vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
}



void plot_two_beams_by_given_alpha_and_phi(double alpha, double phi, bool two_beams = false, bool plot_distortion = false,
                                           double distortion_ampl = 0, double distortion_k = 2*M_PI/lambda)
{
    std::cout << "Plotting for alpha = " << alpha << "\t phi = " << phi;

    double F = get_F_by_beam_parameters_alpha(alpha, phi, beam_width);
    double p = get_p_by_beam_parameters_alpha(alpha, F); // impact parameter по нижней границе
    double h = p + beam_radius;

    std::cout << "\t Focal length = " << F;
    if (plot_distortion) {std::cout << "\t k = " << distortion_k << "\t ampl = " << distortion_ampl;}
    std::cout << std::endl;

    Position p_focus = {0.0, 0.0, 0.0};

    ParabolicSurface mirror1(
        {0.0, 0.0, -F},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        4.0*F, 4.0*F
    );

    ParabolicSurface mirror2(
        {0.0, 0.0, F},
        {-1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        4.0*F, 4.0*F
    );

    PlaneSurface beam_profile({0.0, 0.0, 30.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0});

    Vector distortion_direction = mirror1.dS_over_dxdy({h, 0.0});
    distortion_direction /= distortion_direction.norm();

    double dist_kx =     distortion_k / sqrt(5);
    double dist_ky = 2 * distortion_k / sqrt(5);
    dist_kx /= distortion_direction[2];

    std::vector<SurfaceDistortionHarmonic::DistortionHarmonic> harmonics;

    SurfaceDistortionHarmonic::DistortionHarmonic h1 = { .ampl = distortion_ampl, .kx = dist_kx, .ky = dist_ky};
//    SurfaceDistortion::DistortionHarmonic h2 = {distortion_ampl, dist_kx*1.3, dist_ky*0.09};
//    SurfaceDistortion::DistortionHarmonic h3 = {distortion_ampl, dist_kx*3.2, dist_ky*0.7};
    harmonics.push_back(h1);
//    harmonics.push_back(h2);
//    harmonics.push_back(h3);

    SurfaceDistortionHarmonic mirror1dist(mirror1, distortion_direction, harmonics);

    // Параметры для гауссова пучка, той же дисперсии, что и _|¯¯¯|_, и чтобы в нём энергии как в _|¯¯¯|_
    double sigma = beam_radius / 2;
    double gauss_A = 2;

    ParallelBeamAlpha beam1(lambda, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0},
                            [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 - h) + sqr(x2)) / (2 * sqr(sigma)) ); },
                            //[h](double x1, double x2) { return sqr(x1 - h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
                            //[h](double x1, double x2) { return smoothed(sqrt(sqr(x1 - h) + sqr(x2)), beam_radius); },
                            [](double, double) { return 0.0; });

    ParallelBeamAlpha beam2(lambda, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                            [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 + h) + sqr(x2)) / (2 * sqr(sigma)) ); },
                            //[h](double x1, double x2) { return sqr(x1 + h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
                            //[h](double x1, double x2) { return smoothed(sqrt(sqr(x1 + h) + sqr(x2)), beam_radius); },
                            [](double, double) { return 0.0; });

    double mirror_radius = 3 * sigma;
    //double mirror_radius = beam_radius;
    //double mirror_radius = beam_radius + 0.1;

    SurfaceRegion region1;
    region1.x_min = h - mirror_radius;
    region1.x_max = h + mirror_radius;
    region1.y_min = 0.0 - mirror_radius;
    region1.y_max = 0.0 + mirror_radius;

    SurfaceRegion region2; //           Почему не минус h ???
    region2.x_min = h - mirror_radius;
    region2.x_max = h + mirror_radius;
    region2.y_min = 0.0 - mirror_radius;
    region2.y_max = 0.0 + mirror_radius;

    SurfaceRegion region_profile;
    region_profile.x_min = h - 1.1 * mirror_radius;
    region_profile.x_max = h + 1.1 * mirror_radius;
    region_profile.y_min = 0.0 - 1.1 * mirror_radius;
    region_profile.y_max = 0.0 + 1.1 * mirror_radius;

    StrattonChuReflection reflection1(mirror1, beam1, region1);
    StrattonChuReflection reflection1dist(mirror1dist, beam1, region1);
    StrattonChuReflection reflection2(mirror2, beam2, region2);

    PlaneSurface focal_plane( p_focus, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0} );
    SurfaceRegion focal_region;
    focal_region.x_min = 0.0 - 12 * lambda;
    focal_region.x_max = 0.0 + 12 * lambda;
    focal_region.y_min = 0.0 - 12 * lambda;
    focal_region.y_max = 0.0 + 12 * lambda;
    int focal_points = 301;

    Vector direction_vector = p_focus - mirror1.point({h, 0.0});
    direction_vector /= direction_vector.norm();

    PlaneSurface focal_plane_transversal( p_focus, {0.0, 1.0, 0.0}, Vector(0.0, 1.0, 0.0) % direction_vector );
    SurfaceRegion focal_region_transversal;
    focal_region_transversal.x_min = 0.0 - 10 * lambda;
    focal_region_transversal.x_max = 0.0 + 10 * lambda;
    focal_region_transversal.y_min = 0.0 - 10 * lambda;
    focal_region_transversal.y_max = 0.0 + 10 * lambda;
    int focal_points_transversal = 201;


    std::string filename_suffix = "alpha-";
    filename_suffix += std::to_string(int(alpha*100));
    //std::string filename_suffix = "k-";
    //filename_suffix += std::to_string(int(distortion_k*10));
    //std::string filename_suffix = "ampl-";
    //filename_suffix += std::to_string(int(distortion_ampl*1.0e7));


    plot_field_on_given_surface(beam_profile, beam1, region_profile, 300, "E", "beam_profile", filename_suffix);
    std::cout << "beam profile plotted" << std::endl;

    plot_field_on_given_surface(mirror1, beam1, region1, 100, "E", "mirror1", filename_suffix);
    std::cout << "mirror1 plotted" << std::endl;
    if (plot_distortion)
    {
        plot_field_on_given_surface(mirror1dist, beam1, region1, 100, "E", "mirror1_dist", filename_suffix);
        std::cout << "mirror1dist plotted" << std::endl;
    }
    if (two_beams)
    {
        plot_field_on_given_surface(mirror2, beam2, region2, 100, "E", "mirror2", filename_suffix);
        std::cout << "mirror2 plotted" << std::endl;
    }

    plot_field_on_given_surface(focal_plane, reflection1, focal_region, focal_points, "E1", "longitudinal_E1", filename_suffix);
    std::cout << "longitudinal_E1 plotted" << std::endl;
    if (plot_distortion)
    {
        plot_field_on_given_surface(focal_plane, reflection1dist, focal_region, focal_points, "E1", "longitudinal_E1_dist", filename_suffix);
        std::cout << "longitudinal_E1_dist plotted" << std::endl;
    }
    if (two_beams)
    {
        plot_field_on_given_surface(focal_plane, reflection2, focal_region, focal_points, "E2", "longitudinal_E2", filename_suffix);
        std::cout << "longitudinal_E2 plotted" << std::endl;
    }

    plot_field_on_given_surface(focal_plane_transversal, reflection1, focal_region_transversal, focal_points_transversal, "E1", "transversal_E1", filename_suffix);
    std::cout << "transversal_E1 plotted" << std::endl;
    if (plot_distortion)
    {
        plot_field_on_given_surface(focal_plane_transversal, reflection1dist, focal_region_transversal, focal_points_transversal, "E1", "transversal_E1_dist", filename_suffix);
        std::cout << "transversal_E1_dist plotted" << std::endl;
    }
    if (two_beams)
    {
        plot_field_on_given_surface(focal_plane_transversal, reflection2, focal_region_transversal, focal_points_transversal, "E2", "transversal_E2", filename_suffix);
        std::cout << "transversal_E2 plotted" << std::endl;
    }
}


int main()
{

//     //Код для зависимости картин от альфа
    
//     //const double F = 20.0; // cm
//     const double phi = M_PI / 3;

//     std::cout << "Running stratton-chu computation" << std::endl;
//     int steps_count = 50;

// //    double ampl = 0.1 * lambda;
// //    double k = 2 * M_PI / (beam_radius / 10);
// //    double k = 2 * M_PI / (100 * lambda);

//     double alpha_min = 0.0;
//     double alpha_max = M_PI - phi;
//     //double alpha = M_PI/6;

//     //plot_two_beams_by_given_alpha_and_phi(alpha, phi, 0, 0);


//     for(int i = 40; i < steps_count; i++)
//     {
//         double alpha = alpha_min + (alpha_max - alpha_min) / (steps_count-1) * i;

//         try
//         {
//             plot_two_beams_by_given_alpha_and_phi(alpha, phi, 0, 0);
//         }
//         catch(const std::domain_error& ex)
//         {
//             std::cout << "Exception for alpha = " << alpha << ": " << ex.what() << std::endl;
//         }
//     }
//     //конец этой части кода



    //немонохроматический пучок
    const double phi = M_PI / 3;
    double alpha = M_PI / 3;
    const double omega0 = 2 * M_PI / lambda*c;   
    
    double tauint; //для интенсивности
    tauint=11*pow(10, -15);
    //std::cout << tauint << std::endl;

    double tauf=tauint/(2*sqrt(2*log(2))); //для поля
    //std::cout << tauf << std::endl;

    //funct=exp(-sqr(t) / (4 * sqr(tauf)) * exp(-Complex(0.0, 1.0) * omega0 * t);

    fftw_complex *in, *out;
    fftw_plan plan;
    double T=100*tauf;
    int numb=5000;
    //int numb=625;
    //int numb=100;
    double deltat=T/(numb-1);
    double tstart=-T/10;
    //reinterpret_cast //для совместимости форматов комплексных чисел
    std::complex<double> funct;


    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*numb);
    for (int i = 0; i < numb; i++)
    {
        funct=exp(-sqr(tstart+deltat*i) / (4 * sqr(tauf))) * exp(std::complex<double>(0.0, 1.0) * omega0 * (tstart+deltat*i));
        //funct=exp(-sqr(tstart+deltat*i) / (4 * sqr(tauf))) * cos(omega0 * (tstart+deltat*i));
        in[i][0]=funct.real();
        in[i][1]=funct.imag();
    }
    /*
    //тест
    std::ofstream inr("input.txt");
	for (int i = 0; i < numb; i++) {
        if (in[i][1]<0)
            inr << in[i][0] <<"-"<<abs(in[i][1])<<"i"<< '\n';
        else
            inr << in[i][0] <<"+"<<in[i][1]<<"i"<< '\n';
	}
	inr.close();

    std::ofstream iuto("inputabs.txt");
	for (int i = 0; i < numb; i++) {
         iuto << sqrt(sqr(in[i][0])+sqr(in[i][1])) << '\n';
	}
	iuto.close();
    */
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*numb);
    plan = fftw_plan_dft_1d (numb, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute (plan);

    std::complex<double> *outc=reinterpret_cast<std::complex<double>*>(out);

    // //тест
    // std::ofstream outr("output.txt");
	// for (int i = 0; i < numb; i++) {
    //     if (out[i][1]<0)
    //         outr << out[i][0] <<"-"<<abs(out[i][1])<<"i"<< '\n';
    //     else
    //         outr << out[i][0] <<"+"<<out[i][1]<<"i"<< '\n';
	// }
	// outr.close();

    // std::ofstream outre("output1.txt");
	// for (int i = 0; i < numb; i++) {
    //     outre << outc[i] << '\n';
	// }
	// outre.close();

    std::ofstream outo("outputabs.txt");
	for (int i = 0; i < numb; i++) {
         //outo << sqrt(sqr(out[i][0])+sqr(out[i][1])) << '\n';
         outo << abs(outc[i]) << '\n';
         //std::cout << sqrt(sqr(out[i][0])+sqr(out[i][1])) << std::endl;
	}
	outo.close();



    int central_harmonic=round(T*omega0/2/M_PI);
    double halfshir=M_PI/tauf;
    int left_harmonics_border=round(T*(omega0-halfshir)/2/M_PI);
    int right_harmonics_border=round(T*(omega0+halfshir)/2/M_PI);

    //исправления для учёта всех гармоник
    int harmonics_number_full=right_harmonics_border-left_harmonics_border;
    std::cout << harmonics_number_full << std::endl;
    if (harmonics_number_full==harmonics_number) {
        std::cout << "Calculations are made for all points in spectrum" << std::endl;
    } else {
        std::cout << "Set harmonics_number in constants.hpp to " << harmonics_number_full << std::endl;
    }

    //для учёта избранного числа гармоник
    int harmonic_interval=round(static_cast<double>(right_harmonics_border-left_harmonics_border)/(harmonics_number+1));
    int middle_index=harmonics_number/2;




    //std::cout << central_harmonic << std::endl;


       // Параметры для гауссова пучка, той же дисперсии, что и _|¯¯¯|_, и чтобы в нём энергии как в _|¯¯¯|_
    double sigma = beam_radius / 2;
    //double gauss_A = 2/sqrt(harmonics_number);
    double gauss_A = 2;
    //std::vector<double> gauss_AA(harmonics_number);

    //для учёта избранного числа гармоник
    std::vector<double> freqs(harmonics_number);
    std::vector<double> lambdas(harmonics_number);
    std::vector<std::complex<double>> ampls(harmonics_number);
    freqs[middle_index]=omega0;
    lambdas[middle_index]=lambda;
    ampls[middle_index]=outc[central_harmonic];

    //gauss_AA[middle_index]=2/sqrt(harmonics_number);

    //исправления для учёта всех гармоник
    // std::vector<double> freqs(harmonics_number_full);
    // std::vector<double> lambdas(harmonics_number_full); 
    // std::vector<std::complex<double>> ampls(harmonics_number_full);
    //std::vector<double> freqs_full(harmonics_number_full);
    //std::vector<double> lambdas_full(harmonics_number_full); 
    //std::vector<std::complex<double>> ampls_full(harmonics_number_full);

    // for (int i = 0; i < harmonics_number_full; i++) {
    //     freqs[i]=(left_harmonics_border+i)/T*2*M_PI;
    //     lambdas[i]=2*M_PI*c/freqs[i];
    //     ampls[i]=outc[left_harmonics_border+i];
    //     //std::cout << freqs[i] << std::endl;
    // }

    //freqs[central_harmonic-left_harmonics_border]=omega0;
    //lambdas[central_harmonic-left_harmonics_border]=lambda;


    if (middle_index!=0) {
        for (int i = 1; i <= middle_index; i++) {
            freqs[middle_index+i]=(central_harmonic+harmonic_interval*i)/T*2*M_PI;
            lambdas[middle_index+i]=2*M_PI*c/freqs[middle_index+i];
            ampls[middle_index+i]=outc[central_harmonic+harmonic_interval*i];
            freqs[middle_index-i]=(central_harmonic-harmonic_interval*i)/T*2*M_PI;
            lambdas[middle_index-i]=2*M_PI*c/freqs[middle_index-i];
            ampls[middle_index-i]=outc[central_harmonic-harmonic_interval*i];
            
            //std::cout << freqs[i] << std::endl;

            //gauss_AA[middle_index+i]=2/sqrt(harmonics_number)*abs(outc[central_harmonic])/abs(outc[central_harmonic+harmonic_interval*i]);
            //gauss_AA[middle_index-i]=2/sqrt(harmonics_number)*abs(outc[central_harmonic])/abs(outc[central_harmonic-harmonic_interval*i]);
        }
    }



    // for (int i = 0; i < harmonics_number; i++) {
    //         std::cout << gauss_AA[i] << std::endl;
    // }
    
    

    double F = get_F_by_beam_parameters_alpha(alpha, phi, beam_width);
    double p = get_p_by_beam_parameters_alpha(alpha, F); // impact parameter по нижней границе
    double h = p + beam_radius;

    std::cout << "\t Focal length = " << F<< std::endl;
    /*if (plot_distortion) {std::cout << "\t k = " << distortion_k << "\t ampl = " << distortion_ampl;}
    std::cout << std::endl;*/

    Position p_focus = {0.0, 0.0, 0.0};

    ParabolicSurface mirror1(
        {0.0, 0.0, -F},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        4.0*F, 4.0*F
    );

    ParabolicSurface mirror2(
        {0.0, 0.0, F},
        {-1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        4.0*F, 4.0*F
    );

    PlaneSurface beam_profile({0.0, 0.0, 30.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0});


    //std::cout << gauss_A << std::endl;

    double mirror_radius = 3 * sigma;
    //double mirror_radius = beam_radius;
    //double mirror_radius = beam_radius + 0.1;

    SurfaceRegion region1;
    region1.x_min = h - mirror_radius;
    region1.x_max = h + mirror_radius;
    region1.y_min = 0.0 - mirror_radius;
    region1.y_max = 0.0 + mirror_radius;

    SurfaceRegion region2; //           Почему не минус h ???
    region2.x_min = h - mirror_radius;
    region2.x_max = h + mirror_radius;
    region2.y_min = 0.0 - mirror_radius;
    region2.y_max = 0.0 + mirror_radius;

    SurfaceRegion region_profile;
    region_profile.x_min = h - 1.1 * mirror_radius;
    region_profile.x_max = h + 1.1 * mirror_radius;
    region_profile.y_min = 0.0 - 1.1 * mirror_radius;
    region_profile.y_max = 0.0 + 1.1 * mirror_radius;

    std::vector<ParallelBeamAlpha> beams;
    std::vector<StrattonChuReflection> refs;
    beams.reserve(harmonics_number);
    refs.reserve(harmonics_number);
    for (size_t i=0;i<harmonics_number;i++)
    {
        //double gauss_A=gauss_AA[i];
        //std::cout << gauss_A << "  "<< gauss_AA[i]<< std::endl;
         beams.emplace_back(ParallelBeamAlpha(lambdas[i], {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0},
                            [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 - h) + sqr(x2)) / (2 * sqr(sigma)) ); },
                            //[h](double x1, double x2) { return sqr(x1 - h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
                            //[h](double x1, double x2) { return smoothed(sqrt(sqr(x1 - h) + sqr(x2)), beam_radius); },
                            [](double, double) { return 0.0; }));
         refs.emplace_back(mirror1, beams[i], region1);
    }

    /*
    std::cout << refs.size() << std::endl;
    std::cout << (sizeof(freqs)/sizeof(*freqs)) << std::endl;
    std::cout << (sizeof(ampls)/sizeof(*ampls)) << std::endl;
    */


    SpecIFT impulse(refs, freqs, ampls);

     
    PlaneSurface focal_plane( p_focus, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0} );
    SurfaceRegion focal_region;
    focal_region.x_min = 0.0 - 12 * lambda;
    focal_region.x_max = 0.0 + 12 * lambda;
    focal_region.y_min = 0.0 - 12 * lambda;
    focal_region.y_max = 0.0 + 12 * lambda;


    Vector direction_vector = p_focus - mirror1.point({h, 0.0});
    direction_vector /= direction_vector.norm();

    PlaneSurface focal_plane_transversal( p_focus, {0.0, 1.0, 0.0}, Vector(0.0, 1.0, 0.0) % direction_vector );
    SurfaceRegion focal_region_transversal;
    focal_region_transversal.x_min = 0.0 - 10 * lambda;
    focal_region_transversal.x_max = 0.0 + 10 * lambda;
    focal_region_transversal.y_min = 0.0 - 10 * lambda;
    focal_region_transversal.y_max = 0.0 + 10 * lambda;


    double Tmax=5*tauf;
    int number_of_time_points=25;
    double time_interval=Tmax/number_of_time_points;


    std::cout << harmonics_number << " harmonics modeling starts" << std::endl;
    
    for (size_t i=20;i<number_of_time_points;i++)
    {   
        std::string filename_suffix = "time-";
        filename_suffix += std::to_string(i);
        std::cout <<"file" << i  << std::endl;
        //std::cout << "freqs "<< (sizeof(freqs)/sizeof(*freqs)) << std::endl;
        //std::cout << "ampls "<< (sizeof(ampls)/sizeof(*ampls)) << std::endl;
        plot_field_on_given_surface_with_time(focal_plane, impulse, focal_region, focal_points, i*time_interval, "E1", "longitudinal_E1", filename_suffix);
        std::cout << "longitudinal_E1 plotted" << std::endl;
    
        plot_field_on_given_surface_with_time(focal_plane_transversal, impulse, focal_region_transversal, focal_points_transversal, i*time_interval, "E1", "transversal_E1", filename_suffix);
        std::cout << "transversal_E1 plotted" << std::endl;

    }

    
    std::cout << harmonics_number << " harmonics modeling is over" << std::endl;

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
  
    return 0;

}
