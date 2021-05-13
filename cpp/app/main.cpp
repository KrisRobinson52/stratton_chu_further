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
        const std::string& quantity_name, const std::string& filaname_prefix, const std::string& filename_suffix)
{
    VTKSurfaceSaver vtk_saver(n_points, n_points, quantity_name.c_str());
    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
            Position p = surface.point(xy);
            FieldValue field_value = field.get(p);
            vtk_saver.set_point(i, j, p, vec_modulus(field_value.E));
        }
    }
    vtk_saver.save((filaname_prefix + filename_suffix).c_str());
}

void plot_focal_fields(double x_from, double x_to, double z_from, double z_to, int n_points,
                       IField& field1, IField& field2, std::string filename_suffix)
{
    VTKSurfaceSaver vtk_saver_1r(n_points, n_points, "E1_real");
    VTKSurfaceSaver vtk_saver_1i(n_points, n_points, "E1_imag");
    VTKSurfaceSaver vtk_saver_1n(n_points, n_points, "E1_norm");
    VTKSurfaceSaver vtk_saver_2r(n_points, n_points, "E2_real");
    VTKSurfaceSaver vtk_saver_2i(n_points, n_points, "E2_imag");
    VTKSurfaceSaver vtk_saver_2n(n_points, n_points, "E2_norm");
    VTKSurfaceSaver vtk_saver_r(n_points, n_points, "E_real");
    VTKSurfaceSaver vtk_saver_i(n_points, n_points, "E_imag");
    VTKSurfaceSaver vtk_saver_n(n_points, n_points, "E_norm");

    tbb::parallel_for( int(0), n_points,
        [x_from, x_to, z_from, z_to, n_points, &field1, &field2,
         &vtk_saver_1r, &vtk_saver_1i, &vtk_saver_1n,
         &vtk_saver_2r, &vtk_saver_2i, &vtk_saver_2n,
         &vtk_saver_r,  &vtk_saver_i,  &vtk_saver_n]
        (int i) {
            for (int j = 0; j < n_points; j++)
            {
                Position p(x_from + (x_to - x_from) / (n_points-1) * i, 0.0, z_from + (z_to - z_from) / (n_points-1) * j);
                FieldValue val1 = field1.get(p);
                FieldValue val2 = field2.get(p);
                FieldValue val = val1 + val2;

                std::cout << "Point " << p.str() << ": E = " << val.E.str() << std::endl;

                vtk_saver_1r.set_point(i, j, p, vec_real(val1.E));
                vtk_saver_1i.set_point(i, j, p, vec_imag(val1.E));
                vtk_saver_1n.set_point(i, j, p, vec_modulus(val1.E));

                vtk_saver_2r.set_point(i, j, p, vec_real(val2.E));
                vtk_saver_2i.set_point(i, j, p, vec_imag(val2.E));
                vtk_saver_2n.set_point(i, j, p, vec_modulus(val2.E));

                vtk_saver_r.set_point(i, j, p, vec_real(val.E));
                vtk_saver_i.set_point(i, j, p, vec_imag(val.E));
                vtk_saver_n.set_point(i, j, p, vec_modulus(val.E));
            }
        }
    );

    vtk_saver_1r.save(("field-E1-real-" + filename_suffix).c_str());
    vtk_saver_1i.save(("field-E1-imag-" + filename_suffix).c_str());
    vtk_saver_1n.save(("field-E1-norm-" + filename_suffix).c_str());
    vtk_saver_2r.save(("field-E2-real-" + filename_suffix).c_str());
    vtk_saver_2i.save(("field-E2-imag-" + filename_suffix).c_str());
    vtk_saver_2n.save(("field-E2-norm-" + filename_suffix).c_str());
    vtk_saver_r.save(("field-E-real-" + filename_suffix).c_str());
    vtk_saver_i.save(("field-E-imag-" + filename_suffix).c_str());
    vtk_saver_n.save(("field-E-norm-" + filename_suffix).c_str());
}

void plot_two_beams_by_given_alpha_and_phi(double alpha, double phi)
{
    std::cout << "Plotting for alpha = " << alpha;

    double F = get_F_by_beam_parameters_alpha(alpha, phi, beam_width);
    double p = get_p_by_beam_parameters_alpha(alpha, F); // impact parameter
    double h = p + beam_radius;

    std::cout << "\tFocal length = " << F << std::endl;

    ParabolicSurface surface1(
        {0.0, 0.0, -F},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        4.0*F, 4.0*F
    );

    ParabolicSurface surface2(
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

    StrattonChuReflection reflection1(surface1, beam1, region1);
    StrattonChuReflection reflection2(surface2, beam2, region2);

    int n_points = 301;
    double z_from = 0.0 - 15 * lambda;
    double z_to = 0.0 + 15 * lambda;
    double x_from = 0.0 - 15 * lambda;
    double x_to = 0.0 + 15 * lambda;

    std::string filename_suffix = "alpha-";
    filename_suffix += std::to_string(int(alpha*100));

    plot_focal_fields(x_from, x_to, z_from, z_to, n_points, reflection1, reflection2, filename_suffix);

    SurfaceRegion reflector_region;
    reflector_region.x_min = - h - beam_radius;
    reflector_region.x_max = h + beam_radius;
    reflector_region.y_min = 0.0 - beam_radius;
    reflector_region.y_max = 0.0 + beam_radius;

    plot_field_on_given_surface(surface1, beam1 + beam2, reflector_region, 100, "incident_field", "reflector-", filename_suffix);
}

void plot_two_beams_by_given_h_and_F(double h, double F)
{
    std::cout << "Plotting for h = " << h << std::endl;

    ParabolicSurface surface(
        {0.0, 0.0, -F},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        4.0*F, 4.0*F
    );

    ParallelBeamZ beam1(
        lambda,
        //[h, beam_radius](double x, double y) { return exp( - (sqr(x - h) + sqr(y)) / sqr(beam_radius/2) / 2); },
        [h](double x, double y) { return sqr(x - h) + sqr(y) <= sqr(beam_radius) ? 1.0 : 0.0; },
        [](double, double) { return 0.0; },
        0.0
    );

    ParallelBeamZ beam2(
        lambda,
        //[h, beam_radius](double x, double y) { return exp( - (sqr(x - h) + sqr(y)) / sqr(beam_radius/2) / 2); },
        [h](double x, double y) { return sqr(x + h) + sqr(y) <= sqr(beam_radius) ? -1.0 : 0.0; },
        [](double, double) { return 0.0; },
        0.0
    );

    SurfaceRegion region1;
    region1.x_min = h - beam_radius;
    region1.x_max = h + beam_radius;
    region1.y_min = 0.0 - beam_radius;
    region1.y_max = 0.0 + beam_radius;

    SurfaceRegion region2;
    region2.x_min = - h - beam_radius;
    region2.x_max = - h + beam_radius;
    region2.y_min = 0.0 - beam_radius;
    region2.y_max = 0.0 + beam_radius;

    StrattonChuReflection reflection1(surface, beam1, region1);
    StrattonChuReflection reflection2(surface, beam2, region2);

    int n_points = 201;
    double z_from = 0.0 - 15 * lambda;
    double z_to = 0.0 + 15 * lambda;
    double x_from = 0.0 - 15 * lambda;
    double x_to = 0.0 + 15 * lambda;

    std::string filename_suffix = "h-";
    filename_suffix += std::to_string(int(h*10));

    plot_focal_fields(x_from, x_to, z_from, z_to, n_points, reflection1, reflection2, filename_suffix);
}

int main()
{
    //tbb::task_scheduler_init init(1);

    //const double F = 20.0; // cm
    const double phi = M_PI / 3;

    std::cout << "Running stratton-chu computation" << std::endl;
    //int steps_count = 20;

    double alpha_min = 0.0;
    double alpha_max = M_PI - phi;

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
