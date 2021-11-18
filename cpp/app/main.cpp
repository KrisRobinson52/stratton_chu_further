#include "stratton-chu/plane-surface.hpp"
#include "stratton-chu/parabolic-surface.hpp"
#include "stratton-chu/distorted-surface.hpp"
#include "stratton-chu/parallel-beam.hpp"
#include "stratton-chu/stratton-chu-field.hpp"
#include "stratton-chu/utils.hpp"

#include "stratton-chu/csv-saver.hpp"
#include "stratton-chu/vtk-saver.hpp"

#include <iostream>
#include "thread_pool.hpp"


const double lambda = 0.0000091; // 910 nm
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


    TaskQueue tasks;
    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < n_points; j++)
        {
            tasks.push([&, i, j](){
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
    threads.run();

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
    //const double F = 20.0; // cm
    const double phi = M_PI / 3;

    std::cout << "Running stratton-chu computation" << std::endl;
    int steps_count = 50;

//    double ampl = 0.1 * lambda;
//    double k = 2 * M_PI / (beam_radius / 10);
//    double k = 2 * M_PI / (100 * lambda);

    double alpha_min = 0.0;
    double alpha_max = M_PI - phi;
    //double alpha = M_PI/6;

    //plot_two_beams_by_given_alpha_and_phi(alpha, phi, 0, 0);


    for(int i = 0; i < steps_count; i++)
    {
        double alpha = alpha_min + (alpha_max - alpha_min) / (steps_count-1) * i;

        try
        {
            plot_two_beams_by_given_alpha_and_phi(alpha, phi, 0, 0);
        }
        catch(const std::domain_error& ex)
        {
            std::cout << "Exception for alpha = " << alpha << ": " << ex.what() << std::endl;
        }
    }

    return 0;
}
