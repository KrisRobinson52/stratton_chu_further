#include "stratton-chu/plane-surface.hpp"
#include "stratton-chu/parabolic-surface.hpp"
#include "stratton-chu/parallel-beam.hpp"
#include "stratton-chu/stratton-chu-field.hpp"
#include "stratton-chu/utils.hpp"

#include "stratton-chu/csv-saver.hpp"
#include "stratton-chu/vtk-saver.hpp"

#include <tbb/tbb.h>
#include <iostream>

const double lambda = 0.000091; // 910 nm
const double beam_radius = 10; // cm
const double beam_width = 2 * beam_radius; // cm

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

    tbb::parallel_for( int(0), n_points,
        [&region, n_points, &field, &surface,
         &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m]
        (int i) {
            for (int j = 0; j < n_points; j++)
            {
                Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                Position p = surface.point(xy);
                FieldValue field_value = field.get(p);

                vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
                vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
                vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

                FieldValue field_rotated;
                field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
                field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

                vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
                vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
            }
        }
    );

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

    vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
}



void plot_two_beams_by_given_alpha_and_phi(double alpha, double phi)
{
    std::cout << "Plotting for alpha = " << alpha;

    double F = get_F_by_beam_parameters_alpha(alpha, phi, beam_width);
    double p = get_p_by_beam_parameters_alpha(alpha, F); // impact parameter
    double h = p + beam_radius;

    std::cout << "\tFocal length = " << F << std::endl;

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

    ParallelBeamAlpha beam1(lambda, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0},
                            [h](double x1, double x2) { return exp( - (sqr(x1 - h) + sqr(x2)) / (2 * sqr(beam_radius/4)) ); },
                            //[h](double x1, double x2) { return sqr(x1 - h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
                            [](double, double) { return 0.0; });

    ParallelBeamAlpha beam2(lambda, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                            [h](double x1, double x2) { return exp( - (sqr(x1 + h) + sqr(x2)) / (2 * sqr(beam_radius/4)) ); },
                            //[h](double x1, double x2) { return sqr(x1 + h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
                            [](double, double) { return 0.0; });

    SurfaceRegion region1;
    region1.x_min = h - beam_radius;
    region1.x_max = h + beam_radius;
    region1.y_min = 0.0 - beam_radius;
    region1.y_max = 0.0 + beam_radius;

    SurfaceRegion region2;
    region2.x_min = h - beam_radius;
    region2.x_max = h + beam_radius;
    region2.y_min = 0.0 - beam_radius;
    region2.y_max = 0.0 + beam_radius;

    StrattonChuReflection reflection1(mirror1, beam1, region1);
    StrattonChuReflection reflection2(mirror2, beam2, region2);

    PlaneSurface focal_plane( p_focus, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0} );
    SurfaceRegion focal_region;
    focal_region.x_min = 0.0 - 15 * lambda;
    focal_region.x_max = 0.0 + 15 * lambda;
    focal_region.y_min = 0.0 - 15 * lambda;
    focal_region.y_max = 0.0 + 15 * lambda;
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

    plot_field_on_given_surface(mirror1, beam1 + beam2, region1, 100, "E", "mirror1", filename_suffix);
    std::cout << "mirror1 plotted" << std::endl;
    plot_field_on_given_surface(mirror2, beam1 + beam2, region2, 100, "E", "mirror2", filename_suffix);
    std::cout << "mirror2 plotted" << std::endl;

    plot_field_on_given_surface(focal_plane, reflection1, focal_region, focal_points, "E1", "longitudinal_E1", filename_suffix);
    std::cout << "longitudinal_E1 plotted" << std::endl;
    //plot_field_on_given_surface(focal_plane, reflection2, focal_region, focal_points, "E2", "longitudinal_E2", filename_suffix);
    //std::cout << "longitudinal_E2 plotted" << std::endl;

    plot_field_on_given_surface(focal_plane_transversal, reflection1, focal_region_transversal, focal_points_transversal, "E1", "transversal_E1", filename_suffix);
    std::cout << "transversal_E1 plotted" << std::endl;
    //plot_field_on_given_surface(focal_plane_transversal, reflection2, focal_region_transversal, focal_points_transversal, "E2", "transversal_E2", filename_suffix);
    //std::cout << "transversal_E2 plotted" << std::endl;

}


int main()
{
    tbb::task_scheduler_init init(8);

    //const double F = 20.0; // cm
    const double phi = M_PI / 3;

    std::cout << "Running stratton-chu computation" << std::endl;
    //int steps_count = 20;

    double alpha_min = 0.0;
    double alpha_max = M_PI - phi;

    //plot_two_beams_by_given_alpha_and_phi(M_PI_4, phi);

    for(double alpha = alpha_min; alpha < alpha_max; alpha += M_PI / 10)
    {
        try
        {
            plot_two_beams_by_given_alpha_and_phi(alpha, phi);
        }
        catch(const std::domain_error& ex)
        {
            std::cout << "Exception for alpha = " << alpha << ": " << ex.what() << std::endl;
        }
    }

    return 0;
}
