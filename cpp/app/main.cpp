#include "stratton-chu/plane-surface.hpp"
#include "stratton-chu/parabolic-surface.hpp"
#include "stratton-chu/parallelepiped-volume.hpp"
#include "stratton-chu/distorted-surface.hpp"
#include "stratton-chu/parallel-beam.hpp"
#include "stratton-chu/stratton-chu-field.hpp"
#include "stratton-chu/utils.hpp"
#include "stratton-chu/spec-inverse-fourier.hpp"
#include "stratton-chu/fields-in-region.hpp"
#include "stratton-chu/fields-in-3dregion.hpp"

#include "stratton-chu/csv-saver.hpp"
// #include "stratton-chu/vtk-saver.hpp"

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

#define WITH_NAME(var) var, #var
#define GET_VARIABLE_NAME(var) (#var)

const double lambda = 0.000091; // cm = 910 nm
const double beam_radius = 10; // cm
const double beam_width = 2 * beam_radius; // cm

//ffff

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

// void plot_field_on_given_surface(
//         const ISurface& surface,
//         const IField& field,
//         const SurfaceRegion& region,
//         // double phase,
//         int n_points,
//         const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
// {
//     VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
//     VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());

//     // В тау-тау-эн;
//     VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

//     using Threads = ThreadPool<std::function<void()>>;
//     Threads::TaskQueue tasks;

//     for (int i = 0; i < n_points; i++)
//     {
//         for (int j = 0; j < n_points; j++)
//         {
//             tasks.push([&region, &field, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m, n_points, i, j](){
//                 Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
//                 Position p = surface.point(xy);
//                 FieldValue field_value = field.get(p);
//                 //std::cout << "point ready" << std::endl;
//                 // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

//                 vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
//                 vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
//                 vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

//                 FieldValue field_rotated;
//                 field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
//                 field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//                 field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

//                 vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
//                 vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
//             });
//         }
//     }
//     Threads threads{tasks};
//     threads.run();

//     vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
//     vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
//     vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

//     vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
//     vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
// }

// void plot_field_on_given_surface_mod(
//         const ISurface& surface,
//         const IField& field,
//         const SurfaceRegion& region,
//         double phase,
//         int n_points,
//         const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
// {
//     VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
//     VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());

//     // В тау-тау-эн;
//     VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

//     using Threads = ThreadPool<std::function<void()>>;
//     Threads::TaskQueue tasks;

//     for (int i = 0; i < n_points; i++)
//     {
//         for (int j = 0; j < n_points; j++)
//         {
//             tasks.push([&region, &field, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m, n_points, i, j, phase](){
//                 Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
//                 Position p = surface.point(xy);
//                 FieldValue field_value = field.get(p);
//                 //std::cout << "point ready" << std::endl;
//                 field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

//                 vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
//                 vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
//                 vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

//                 FieldValue field_rotated;
//                 field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
//                 field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//                 field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

//                 vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
//                 vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
//             });
//         }
//     }
//     Threads threads{tasks};
//     threads.run();

//     vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
//     vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
//     vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

//     vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
//     vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
// }

// //с учётом времени
// void plot_field_on_given_surface_with_time(
//         const ISurface& surface,
//         const SpecIFT& field,
//         const SurfaceRegion& region,
//         int n_points,
//         double timing,
//         const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
// {
//     VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
//     VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());

//     // В тау-тау-эн;
//     VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

//     using Threads = ThreadPool<std::function<void()>>;
//     Threads::TaskQueue tasks;

//     for (int i = 0; i < n_points; i++)
//     {
//         for (int j = 0; j < n_points; j++)
//         {
//             tasks.push([&region, &field, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m, n_points, timing, i, j](){
//                 Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
//                 Position p = surface.point(xy);
//                 FieldValue field_value = field.get(p, timing);
//                 //std::cout << "point ready "<< i << "," << j << std::endl;

//                 vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
//                 vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
//                 vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

//                 FieldValue field_rotated;
//                 field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
//                 field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//                 field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

//                 vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
//                 vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
//             });
//         }
//     }
//     Threads threads{tasks};
//     threads.run();

//     vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
//     vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
//     vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

//     vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
//     vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
// }


// //с учётом времени (3d array method)
// void plot_field_on_given_surface_with_time_3darr(
//         const ISurface& surface,
//         const SurfaceRegion& region,
//         auto& fields, //FieldsInRegion<>& fields;
//         std::vector<double> &freqs,
//         std::vector<std::complex<double>> &amplF,
//         int n_points,
//         double timing,
//         const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
// {
//     VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
//     VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());

//     // В тау-тау-эн;
//     VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

//     fields.constructSummary(freqs, amplF, timing);

//     using Threads = ThreadPool<std::function<void()>>;
//     Threads::TaskQueue tasks;

//     for (int i = 0; i < n_points; i++)
//     {
//         for (int j = 0; j < n_points; j++)
//         {
//             tasks.push([&region, &fields, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m, n_points, timing, i, j](){
//                 Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
//                 Position p = surface.point(xy);
//                 FieldValue field_value = fields.getSummary(i, j);
//                 //std::cout << "point ready "<< i << "," << j << std::endl;

//                 vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
//                 vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
//                 vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

//                 FieldValue field_rotated;
//                 field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
//                 field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//                 field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

//                 vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
//                 vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
//             });
//         }
//     }
//     Threads threads{tasks};
//     threads.run();

//     vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
//     vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
//     vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

//     vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
//     vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
// }

// //3d array method
// void plot_field_on_given_surface_3darr(
//         const ISurface& surface,
//         const SurfaceRegion& region,
//         auto& fields, //FieldsInRegion<>& fields;
//         // std::vector<double> &phases,
//         // size_t iter,
//         int n_points,
//         const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
// {
//     VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
//     VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());
//     VTKSurfaceSaver vtk_saver_int(n_points, n_points, (quantity_name + "_intens").c_str());

//     // // В тау-тау-эн;
//     // VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
//     // VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

//     using Threads = ThreadPool<std::function<void()>>;
//     Threads::TaskQueue tasks;

//     for (int i = 0; i < n_points; i++)
//     {
//         for (int j = 0; j < n_points; j++)
//         {
//             tasks.push([&region, &fields, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_int, /*&vtk_saver_surf_m,*/ n_points,  i, j/*, iter*/](){
//                 Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
//                 Position p = surface.point(xy);
//                 FieldValue field_value = fields.getSummary(i, j);
//                 // FieldValue field_value = fields.getHarmonic(iter, i, j);
//                 //std::cout << "point ready "<< i << "," << j << std::endl;

//                 vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
//                 vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
//                 vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

//                 double absE;
//                 Vector ves;
//                 absE=0.0;
//                 ves=max_field(field_value.E);
//                 for (size_t l=0;l<3;l++){
//                     absE+=ves[l]*ves[l];
//                 }
//                 absE=absE*4e13;
//                 ves[0]=absE;
//                 ves[1]=0;
//                 ves[2]=0;
//                 vtk_saver_int.set_point(i, j, p, ves);

//                 // FieldValue field_rotated;
//                 // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
//                 // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//                 // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

//                 // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
//                 // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
//             });
//         }
//     }
//     Threads threads{tasks};
//     threads.run();

//     vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
//     vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
//     vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());
//     vtk_saver_int.save((filename_prefix + "_intens_" + filename_suffix).c_str());

//     // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
//     // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
// }

// void plot_field_on_given_surface_3darr_1(
//         const ISurface& surface,
//         const SurfaceRegion& region,
//         auto& fields, //FieldsInRegion<>& fields;
//         // std::vector<double> &phases,
//         // size_t iter,
//         int n_points,
//         const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
// {
//     VTKSurfaceSaver vtk_saver_r(55, 55, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_i(55, 55, (quantity_name + "_imag").c_str());
//     VTKSurfaceSaver vtk_saver_m(55, 55, (quantity_name + "_max").c_str());
//     VTKSurfaceSaver vtk_saver_int(55, 55, (quantity_name + "_intens").c_str());

//     // // В тау-тау-эн;
//     // VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
//     // VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

//     using Threads = ThreadPool<std::function<void()>>;
//     Threads::TaskQueue tasks;

//     for (int i = 123; i < 178; i++)
//     {
//         for (int j = 123; j < 178; j++)
//         {
//             tasks.push([&region, &fields, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_int, /*&vtk_saver_surf_m,*/ n_points,  i, j/*, iter*/](){
//                 Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
//                 Position p = surface.point(xy);
//                 FieldValue field_value = fields.getSummary(i, j);
//                 // FieldValue field_value = fields.getHarmonic(iter, i, j);
//                 //std::cout << "point ready "<< i << "," << j << std::endl;

//                 vtk_saver_r.set_point(i-123, j-123, p, vec_real(field_value.E));
//                 vtk_saver_i.set_point(i-123, j-123, p, vec_imag(field_value.E));
//                 vtk_saver_m.set_point(i-123, j-123, p, max_field(field_value.E));

//                 double absE;
//                 Vector ves;
//                 absE=0.0;
//                 ves=max_field(field_value.E);
//                 for (size_t l=0;l<3;l++){
//                     absE+=ves[l]*ves[l];
//                 }
//                 absE=absE*4e13*50/16;
//                 ves[0]=absE;
//                 ves[1]=0;
//                 ves[2]=0;
//                 vtk_saver_int.set_point(i-123, j-123, p, ves);

//                 // FieldValue field_rotated;
//                 // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
//                 // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//                 // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

//                 // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
//                 // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
//             });
//         }
//     }
//     Threads threads{tasks};
//     threads.run();

//     vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
//     vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
//     vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());
//     vtk_saver_int.save((filename_prefix + "_intens_" + filename_suffix).c_str());

//     // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
//     // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
// }

// void plot_field_on_given_surface_3darr_2(
//         const ISurface& surface,
//         const SurfaceRegion& region,
//         auto& fields, //FieldsInRegion<>& fields;
//         // std::vector<double> &phases,
//         // size_t iter,
//         int n_points,
//         const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
// {
//     VTKSurfaceSaver vtk_saver_r(45, 45, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_i(45, 45, (quantity_name + "_imag").c_str());
//     VTKSurfaceSaver vtk_saver_m(45, 45, (quantity_name + "_max").c_str());
//     VTKSurfaceSaver vtk_saver_int(45, 45, (quantity_name + "_intens").c_str());

//     // // В тау-тау-эн;
//     // VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
//     // VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

//     using Threads = ThreadPool<std::function<void()>>;
//     Threads::TaskQueue tasks;

//     for (int i = 78; i < 123; i++)
//     {
//         for (int j = 78; j < 123; j++)
//         {
//             tasks.push([&region, &fields, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_int, /*&vtk_saver_surf_m,*/ n_points,  i, j/*, iter*/](){
//                 Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
//                 Position p = surface.point(xy);
//                 FieldValue field_value = fields.getSummary(i, j);
//                 // FieldValue field_value = fields.getHarmonic(iter, i, j);
//                 //std::cout << "point ready "<< i << "," << j << std::endl;

//                 vtk_saver_r.set_point(i-78, j-78, p, vec_real(field_value.E));
//                 vtk_saver_i.set_point(i-78, j-78, p, vec_imag(field_value.E));
//                 vtk_saver_m.set_point(i-78, j-78, p, max_field(field_value.E));

//                 double absE;
//                 Vector ves;
//                 absE=0.0;
//                 ves=max_field(field_value.E);
//                 for (size_t l=0;l<3;l++){
//                     absE+=ves[l]*ves[l];
//                 }
//                 absE=absE*4e13*50/16;
//                 ves[0]=absE;
//                 ves[1]=0;
//                 ves[2]=0;
//                 vtk_saver_int.set_point(i-78, j-78, p, ves);

//                 // FieldValue field_rotated;
//                 // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
//                 // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//                 // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

//                 // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
//                 // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
//             });
//         }
//     }
//     Threads threads{tasks};
//     threads.run();

//     vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
//     vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
//     vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());
//     vtk_saver_int.save((filename_prefix + "_intens_" + filename_suffix).c_str());

//     // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
//     // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
// }

// void plot_field_in_given_area_3darr(
//         const IVolume& volume,
//         const VolumeRegion& region,
//         auto& fields, //FieldsInRegion<>& fields;
//         // std::vector<double> &phases,
//         // size_t iter,
//         int n_points,
//         const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
// {
//     VTKVolumeSaver vtk_saver_r(n_points, n_points, n_points, (quantity_name + "_real").c_str());
//     VTKVolumeSaver vtk_saver_i(n_points, n_points, n_points, (quantity_name + "_imag").c_str());
//     VTKVolumeSaver vtk_saver_m(n_points, n_points, n_points, (quantity_name + "_max").c_str());
//     VTKVolumeSaver vtk_saver_int(n_points, n_points, n_points, (quantity_name + "_intens").c_str());

//     using Threads = ThreadPool<std::function<void()>>;
//     Threads::TaskQueue tasks;

//     for (int k = 0; k < n_points; k++)
//     {
//         for (int i = 0; i < n_points; i++)
//         {
//             for (int j = 0; j < n_points; j++)
//             {
//                 tasks.push([&region, &fields, &volume, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_int, /*&vtk_saver_surf_m,*/ n_points,  i, j, k/*, iter*/](){
//                     Vector xyz(region.x_min + region.widthx() / (n_points-1) * i, region.y_min + region.widthy() / (n_points-1) * j, region.z_min + region.height() / (n_points-1) * k);
//                     Position p = volume.point(xyz);
//                     FieldValue field_value = fields.getSummary(i, j, k);
//                     // FieldValue field_value = fields.getHarmonic(iter, i, j);
//                     //std::cout << "point ready "<< i << "," << j << std::endl;

//                     vtk_saver_r.set_point(i, j, k, p, vec_real(field_value.E));
//                     vtk_saver_i.set_point(i, j, k, p, vec_imag(field_value.E));
//                     vtk_saver_m.set_point(i, j, k, p, max_field(field_value.E));

//                     double absE;
//                     Vector ves;
//                     absE=0.0;
//                     ves=max_field(field_value.E);
//                     for (size_t l=0;l<3;l++){
//                         absE+=ves[l]*ves[l];
//                     }
//                     absE=absE*4e13*50/16;
//                     ves[0]=absE;
//                     ves[1]=0;
//                     ves[2]=0;
//                     vtk_saver_int.set_point(i, j, k, p, ves);

//                     // FieldValue field_rotated;
//                     // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
//                     // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//                     // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

//                     // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
//                     // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
//                 });
//             }
//         }
//     }
//     Threads threads{tasks};
//     threads.run();

//     vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
//     vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
//     vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());
//     vtk_saver_int.save((filename_prefix + "_intens_" + filename_suffix).c_str());

//     // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
//     // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
// }

// void plot_two_beams_by_given_alpha_and_phi(double alpha, double phi, bool two_beams = false, bool plot_distortion = false,
//                                            double distortion_ampl = 0, double distortion_k = 2*M_PI/lambda)
// {
//     std::cout << "Plotting for alpha = " << alpha << "\t phi = " << phi;

//     double F = get_F_by_beam_parameters_alpha(alpha, phi, beam_width);
//     double p = get_p_by_beam_parameters_alpha(alpha, F); // impact parameter по нижней границе
//     double h = p + beam_radius;

//     std::cout << "\t Focal length = " << F;
//     if (plot_distortion) {std::cout << "\t k = " << distortion_k << "\t ampl = " << distortion_ampl;}
//     std::cout << std::endl;

//     Position p_focus = {0.0, 0.0, 0.0};

//     ParabolicSurface mirror1(
//         {0.0, 0.0, -F},
//         {1.0, 0.0, 0.0},
//         {0.0, 1.0, 0.0},
//         4.0*F, 4.0*F
//     );

//     ParabolicSurface mirror2(
//         {0.0, 0.0, F},
//         {-1.0, 0.0, 0.0},
//         {0.0, 1.0, 0.0},
//         4.0*F, 4.0*F
//     );

//     PlaneSurface beam_profile({0.0, 0.0, 30.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0});

//     Vector distortion_direction = mirror1.dS_over_dxdy({h, 0.0});
//     distortion_direction /= distortion_direction.norm();

//     double dist_kx =     distortion_k / sqrt(5);
//     double dist_ky = 2 * distortion_k / sqrt(5);
//     dist_kx /= distortion_direction[2];

//     std::vector<SurfaceDistortionHarmonic::DistortionHarmonic> harmonics;

//     SurfaceDistortionHarmonic::DistortionHarmonic h1 = { .ampl = distortion_ampl, .kx = dist_kx, .ky = dist_ky};
// //    SurfaceDistortion::DistortionHarmonic h2 = {distortion_ampl, dist_kx*1.3, dist_ky*0.09};
// //    SurfaceDistortion::DistortionHarmonic h3 = {distortion_ampl, dist_kx*3.2, dist_ky*0.7};
//     harmonics.push_back(h1);
// //    harmonics.push_back(h2);
// //    harmonics.push_back(h3);

//     SurfaceDistortionHarmonic mirror1dist(mirror1, distortion_direction, harmonics);

//     // Параметры для гауссова пучка, той же дисперсии, что и _|¯¯¯|_, и чтобы в нём энергии как в _|¯¯¯|_
//     // double sigma = beam_radius / 2;
//     double sigma = beam_radius;
//     double gauss_A = 1;
    
//     ParallelBeamAlpha beam1(lambda, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0},
//                             // [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 - h) + sqr(x2)) / (2 * sqr(sigma)) ); },
//                             [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (2 * pow(sigma, 12)) ); },
//                             //[h](double x1, double x2) { return sqr(x1 - h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
//                             //[h](double x1, double x2) { return smoothed(sqrt(sqr(x1 - h) + sqr(x2)), beam_radius); },
//                             [](double, double) { return 0.0; });

//     ParallelBeamAlpha beam2(lambda, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
//                             [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 + h) + sqr(x2)) / (2 * sqr(sigma)) ); },
//                             //[h](double x1, double x2) { return sqr(x1 + h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
//                             //[h](double x1, double x2) { return smoothed(sqrt(sqr(x1 + h) + sqr(x2)), beam_radius); },
//                             [](double, double) { return 0.0; });

//     double mirror_radius = 3 * sigma;
//     //double mirror_radius = beam_radius;
//     //double mirror_radius = beam_radius + 0.1;

//     SurfaceRegion region1;
//     region1.x_min = h - mirror_radius;
//     region1.x_max = h + mirror_radius;
//     region1.y_min = 0.0 - mirror_radius;
//     region1.y_max = 0.0 + mirror_radius;

//     SurfaceRegion region2; //           Почему не минус h ???
//     region2.x_min = h - mirror_radius;
//     region2.x_max = h + mirror_radius;
//     region2.y_min = 0.0 - mirror_radius;
//     region2.y_max = 0.0 + mirror_radius;

//     SurfaceRegion region_profile;
//     region_profile.x_min = h - 1.1 * mirror_radius;
//     region_profile.x_max = h + 1.1 * mirror_radius;
//     region_profile.y_min = 0.0 - 1.1 * mirror_radius;
//     region_profile.y_max = 0.0 + 1.1 * mirror_radius;

//     StrattonChuReflection reflection1(mirror1, beam1, region1);
//     StrattonChuReflection reflection1dist(mirror1dist, beam1, region1);
//     StrattonChuReflection reflection2(mirror2, beam2, region2);

//     PlaneSurface focal_plane( p_focus, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0} );
//     SurfaceRegion focal_region;
//     focal_region.x_min = 0.0 - 12 * lambda;
//     focal_region.x_max = 0.0 + 12 * lambda;
//     focal_region.y_min = 0.0 - 12 * lambda;
//     focal_region.y_max = 0.0 + 12 * lambda;
//     int focal_points = 301;

//     Vector direction_vector = p_focus - mirror1.point({h, 0.0});
//     direction_vector /= direction_vector.norm();

//     PlaneSurface focal_plane_transversal( p_focus, {0.0, 1.0, 0.0}, Vector(0.0, 1.0, 0.0) % direction_vector );
//     SurfaceRegion focal_region_transversal;
//     focal_region_transversal.x_min = 0.0 - 10 * lambda;
//     focal_region_transversal.x_max = 0.0 + 10 * lambda;
//     focal_region_transversal.y_min = 0.0 - 10 * lambda;
//     focal_region_transversal.y_max = 0.0 + 10 * lambda;
//     int focal_points_transversal = 201;


//     std::string filename_suffix = "alpha-";
//     filename_suffix += std::to_string(int(alpha*100));
//     //std::string filename_suffix = "k-";
//     //filename_suffix += std::to_string(int(distortion_k*10));
//     //std::string filename_suffix = "ampl-";
//     //filename_suffix += std::to_string(int(distortion_ampl*1.0e7));


//     plot_field_on_given_surface(beam_profile, beam1, region_profile, 300, "E", "beam_profile", filename_suffix);
//     std::cout << "beam profile plotted" << std::endl;

//     plot_field_on_given_surface(mirror1, beam1, region1, 100, "E", "mirror1", filename_suffix);
//     std::cout << "mirror1 plotted" << std::endl;
//     if (plot_distortion)
//     {
//         plot_field_on_given_surface(mirror1dist, beam1, region1, 100, "E", "mirror1_dist", filename_suffix);
//         std::cout << "mirror1dist plotted" << std::endl;
//     }
//     if (two_beams)
//     {
//         plot_field_on_given_surface(mirror2, beam2, region2, 100, "E", "mirror2", filename_suffix);
//         std::cout << "mirror2 plotted" << std::endl;
//     }

//     plot_field_on_given_surface(focal_plane, reflection1, focal_region, focal_points, "E1", "longitudinal_E1", filename_suffix);
//     std::cout << "longitudinal_E1 plotted" << std::endl;
//     if (plot_distortion)
//     {
//         plot_field_on_given_surface(focal_plane, reflection1dist, focal_region, focal_points, "E1", "longitudinal_E1_dist", filename_suffix);
//         std::cout << "longitudinal_E1_dist plotted" << std::endl;
//     }
//     if (two_beams)
//     {
//         plot_field_on_given_surface(focal_plane, reflection2, focal_region, focal_points, "E2", "longitudinal_E2", filename_suffix);
//         std::cout << "longitudinal_E2 plotted" << std::endl;
//     }

//     plot_field_on_given_surface(focal_plane_transversal, reflection1, focal_region_transversal, focal_points_transversal, "E1", "transversal_E1", filename_suffix);
//     std::cout << "transversal_E1 plotted" << std::endl;
//     if (plot_distortion)
//     {
//         plot_field_on_given_surface(focal_plane_transversal, reflection1dist, focal_region_transversal, focal_points_transversal, "E1", "transversal_E1_dist", filename_suffix);
//         std::cout << "transversal_E1_dist plotted" << std::endl;
//     }
//     if (two_beams)
//     {
//         plot_field_on_given_surface(focal_plane_transversal, reflection2, focal_region_transversal, focal_points_transversal, "E2", "transversal_E2", filename_suffix);
//         std::cout << "transversal_E2 plotted" << std::endl;
//     }
// }

// void plot_two_beams_by_given_alpha_and_phi_mod(double alpha, double p, size_t beams_number, size_t iter, bool two_beams = false, bool plot_distortion = false,
//                                            double distortion_ampl = 0, double distortion_k = 2*M_PI/lambda)
// {
//     std::cout << "Plotting for alpha = " << alpha;

   

    
//     double F = get_F_by_beam_parameters_alpha_and_p(alpha, p);
//     double h = p + beam_radius;

//     std::cout << "\t Focal length = " << F;
//     if (plot_distortion) {std::cout << "\t k = " << distortion_k << "\t ampl = " << distortion_ampl;}
//     std::cout << std::endl;

//     Position p_focus = {0.0, 0.0, 0.0};

//     ParabolicSurface mirror1(
//         {0.0, 0.0, -F},
//         {1.0*cos(2*M_PI/beams_number*iter), 1.0*sin(2*M_PI/beams_number*iter), 0.0},
//         {-1.0*sin(2*M_PI/beams_number*iter), 1.0*cos(2*M_PI/beams_number*iter), 0.0},
//         4.0*F, 4.0*F
//     );

//     ParabolicSurface mirror2(
//         {0.0, 0.0, F},
//         {-1.0*cos(2*M_PI/beams_number*iter), -1.0*sin(2*M_PI/beams_number*iter), 0.0},
//         {-1.0*sin(2*M_PI/beams_number*iter), 1.0*cos(2*M_PI/beams_number*iter), 0.0},
//         4.0*F, 4.0*F
//     );

//     PlaneSurface beam_profile({0.0, 0.0, 30.0}, {1.0*cos(2*M_PI/beams_number*iter), 1.0*sin(2*M_PI/beams_number*iter), 0.0}, {1.0*sin(2*M_PI/beams_number*iter), -1.0*cos(2*M_PI/beams_number*iter), 0.0});

//     Vector distortion_direction = mirror1.dS_over_dxdy({h, 0.0});
//     distortion_direction /= distortion_direction.norm();

//     double dist_kx =     distortion_k / sqrt(5);
//     double dist_ky = 2 * distortion_k / sqrt(5);
//     dist_kx /= distortion_direction[2];

//     std::vector<SurfaceDistortionHarmonic::DistortionHarmonic> harmonics;

//     SurfaceDistortionHarmonic::DistortionHarmonic h1 = { .ampl = distortion_ampl, .kx = dist_kx, .ky = dist_ky};
// //    SurfaceDistortion::DistortionHarmonic h2 = {distortion_ampl, dist_kx*1.3, dist_ky*0.09};
// //    SurfaceDistortion::DistortionHarmonic h3 = {distortion_ampl, dist_kx*3.2, dist_ky*0.7};
//     harmonics.push_back(h1);
// //    harmonics.push_back(h2);
// //    harmonics.push_back(h3);

//     SurfaceDistortionHarmonic mirror1dist(mirror1, distortion_direction, harmonics);

//     // Параметры для гауссова пучка, той же дисперсии, что и _|¯¯¯|_, и чтобы в нём энергии как в _|¯¯¯|_
//     double sigma = beam_radius;
//     double gauss_A = 2;

//     ParallelBeamAlpha beam1(lambda, {0.0, 0.0, 0.0}, {1.0*cos(2*M_PI/beams_number*iter), 1.0*sin(2*M_PI/beams_number*iter), 0.0}, {1.0*sin(2*M_PI/beams_number*iter), -1.0*cos(2*M_PI/beams_number*iter), 0.0},
//                             //  [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 - h) + sqr(x2)) / (2 * sqr(sigma)) ); },
//                             [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
//                             //[h](double x1, double x2) { return sqr(x1 - h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
//                             // [h](double x1, double x2) { return smoothed(sqrt(sqr(x1 - h) + sqr(x2)), beam_radius); },
//                             [](double, double) { return 0.0; });

//     ParallelBeamAlpha beam2(lambda, {0.0, 0.0, 0.0}, {1.0*cos(2*M_PI/beams_number*iter), 1.0*sin(2*M_PI/beams_number*iter), 0.0}, {-1.0*sin(2*M_PI/beams_number*iter), 1.0*cos(2*M_PI/beams_number*iter), 0.0},
//                             //  [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 + h) + sqr(x2)) / (2 * sqr(sigma)) ); },
//                             [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
//                             //[h](double x1, double x2) { return sqr(x1 + h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
//                             // [h](double x1, double x2) { return smoothed(sqrt(sqr(x1 + h) + sqr(x2)), beam_radius); },
//                             [](double, double) { return 0.0; });

//     double mirror_radius = 3 * sigma;
//     //double mirror_radius = beam_radius;
//     //double mirror_radius = beam_radius + 0.1;

//     SurfaceRegion region1;
//     region1.x_min = h - mirror_radius;
//     region1.x_max = h + mirror_radius;
//     region1.y_min = 0.0 - mirror_radius;
//     region1.y_max = 0.0 + mirror_radius;

//     SurfaceRegion region2; //           Почему не минус h ???
//     region2.x_min = h - mirror_radius;
//     region2.x_max = h + mirror_radius;
//     region2.y_min = 0.0 - mirror_radius;
//     region2.y_max = 0.0 + mirror_radius;

//     SurfaceRegion region_profile;
//     region_profile.x_min = h - 1.1 * mirror_radius;
//     region_profile.x_max = h + 1.1 * mirror_radius;
//     region_profile.y_min = 0.0 - 1.1 * mirror_radius;
//     region_profile.y_max = 0.0 + 1.1 * mirror_radius;

//     // refs1.emplace_back(mirror1, beam1, region1);
//     // refs2.emplace_back(mirror2, beam2, region2);

//     StrattonChuReflection reflection1(mirror1, beam1, region1);
//     StrattonChuReflection reflection1dist(mirror1dist, beam1, region1);
//     StrattonChuReflection reflection2(mirror2, beam2, region2);

//     // PlaneSurface focal_plane( p_focus, {0.0, 0.0, 1.0}, {1.0*cos(M_PI/3*iter), 1.0*sin(M_PI/3*iter), 0.0} );
//     PlaneSurface focal_plane( p_focus, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} );
//     SurfaceRegion focal_region;
//     focal_region.x_min = 0.0 - 12 * lambda;
//     focal_region.x_max = 0.0 + 12 * lambda;
//     focal_region.y_min = 0.0 - 12 * lambda;
//     focal_region.y_max = 0.0 + 12 * lambda;
//     int focal_points = 301;

//     // Vector direction_vector = p_focus - mirror1.point({h, 0.0});
//     // direction_vector /= direction_vector.norm();

//     // PlaneSurface focal_plane_transversal( p_focus, {-1.0*sin(M_PI/3*iter), 1.0*cos(M_PI/3*iter), 0.0}, Vector(-1.0*sin(M_PI/3*iter), 1.0*cos(M_PI/3*iter), 0.0) % direction_vector );
//     PlaneSurface focal_plane_transversal( p_focus, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0} );
//     SurfaceRegion focal_region_transversal;
//     focal_region_transversal.x_min = 0.0 - 10 * lambda;
//     focal_region_transversal.x_max = 0.0 + 10 * lambda;
//     focal_region_transversal.y_min = 0.0 - 10 * lambda;
//     focal_region_transversal.y_max = 0.0 + 10 * lambda;
//     int focal_points_transversal = 201;


//     std::string filename_suffix = "alpha-";
//     filename_suffix += std::to_string(int(alpha*100));
//     filename_suffix += std::to_string(iter);
    
//     //std::string filename_suffix = "k-";
//     //filename_suffix += std::to_string(int(distortion_k*10));
//     //std::string filename_suffix = "ampl-";
//     //filename_suffix += std::to_string(int(distortion_ampl*1.0e7));


//     // plot_field_on_given_surface(beam_profile, beam1, region_profile, 300, "E", "beam_profile", filename_suffix);
//     // std::cout << "beam profile plotted" << std::endl;

//     plot_field_on_given_surface(mirror1, beam1, region1, 100, "E", "mirror1", filename_suffix);
//     std::cout << "mirror1 plotted" << std::endl;
//     if (plot_distortion)
//     {
//         plot_field_on_given_surface(mirror1dist, beam1, region1, 100, "E", "mirror1_dist", filename_suffix);
//         std::cout << "mirror1dist plotted" << std::endl;
//     }
//     if (two_beams)
//     {
//         plot_field_on_given_surface(mirror2, beam2, region2, 100, "E", "mirror2", filename_suffix);
//         std::cout << "mirror2 plotted" << std::endl;
//     }
//     double phase1=0;
//     double phase2=0;
//     plot_field_on_given_surface_mod(focal_plane, reflection1, focal_region, phase1, focal_points, "E1", "longitudinal_E1", filename_suffix);
//     std::cout << "longitudinal_E1 plotted" << std::endl;
//     if (plot_distortion)
//     {
//         plot_field_on_given_surface(focal_plane, reflection1dist, focal_region, focal_points, "E1", "longitudinal_E1_dist", filename_suffix);
//         std::cout << "longitudinal_E1_dist plotted" << std::endl;
//     }
//     if (two_beams)
//     {
//         plot_field_on_given_surface_mod(focal_plane, reflection2, focal_region, phase2, focal_points, "E2", "longitudinal_E2", filename_suffix);
//         std::cout << "longitudinal_E2 plotted" << std::endl;
//     }

//     plot_field_on_given_surface_mod(focal_plane_transversal, reflection1, focal_region_transversal, phase1, focal_points_transversal, "E1", "transversal_E1", filename_suffix);
//     std::cout << "transversal_E1 plotted" << std::endl;
//     if (plot_distortion)
//     {
//         plot_field_on_given_surface(focal_plane_transversal, reflection1dist, focal_region_transversal, focal_points_transversal, "E1", "transversal_E1_dist", filename_suffix);
//         std::cout << "transversal_E1_dist plotted" << std::endl;
//     }
//     if (two_beams)
//     {
//         plot_field_on_given_surface_mod(focal_plane_transversal, reflection2, focal_region_transversal, phase2, focal_points_transversal, "E2", "transversal_E2", filename_suffix);
//         std::cout << "transversal_E2 plotted" << std::endl;
//     }
// }




int main(int argc, const char* argv[])
{
    //Код для монохроматических пучков
    
    //const double F = 20.0; // cm
    // const double phi = M_PI / 3;

    std::cout << "Running stratton-chu computation " << argv[1] << ' ' << argv[2] <<std::endl;
    // int steps_count = 50;

//    double ampl = 0.1 * lambda;
//    double k = 2 * M_PI / (beam_radius / 10);
//    double k = 2 * M_PI / (100 * lambda);

    // double alpha_min = 0.0;
    // double alpha_max = M_PI - phi;
    size_t beams_number = std::atoll(argv[1]);

    //double belt
    // double alpha = M_PI/2;
    // // double p = beam_radius*sqrt(3);
    // double p = beam_radius/tan(M_PI/beams_number);    
    // double F = get_F_by_beam_parameters_alpha_and_p(alpha, p);
    // double h = p + beam_radius; 

    //double belt cfg 2

    // double alpha = 35.6*M_PI/180;
    // double phi = 43*M_PI/180; 
    // double gamma_med = M_PI/6; 

    double betta, gamma_med, delta, alpha, phic, phi, deltab;
    // betta= 2*M_PI/beams_number*0.89;
    deltab=0.1;
    betta= 2*M_PI/beams_number-deltab;
    delta = 0.1;
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

    double F = get_F_by_beam_parameters_alpha(alpha, phi, beam_width);
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
    std::cout << std::endl;

    double sigma = beam_radius;
    double gauss_A = 1;
    // double mirror_radius = 3 * sigma;
    double mirror_radius = 1.5 * sigma;


    // std::vector<ParabolicSurface> mirrors;
    // mirrors.reserve(beams_number);
    std::vector<ParabolicSurface> mirrors1;
    mirrors1.reserve(beams_number);
    std::vector<ParabolicSurface> mirrors2;
    mirrors2.reserve(beams_number);

    

    // std::vector<ParallelBeamAlpha> beams;
    // beams.reserve(beams_number);
    std::vector<ParallelBeamAlpha> beams1;
    beams1.reserve(beams_number);
    std::vector<ParallelBeamAlpha> beams2;
    beams2.reserve(beams_number);

    Position p_focus = {0.0, 0.0, 0.0};
    
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

    double add_angle;
    if (beams_number%2==1){
        add_angle=M_PI/beams_number;
        // add_angle=0;
    }
    else {
        // add_angle=M_PI/beams_number;
        add_angle=0;
    }
    
    Vector ax1m1, ax2m1, axFm1;
    Vector ax1m2, ax2m2, axFm2;

    std::vector<double> mirr_adj;
    mirr_adj.reserve(2*beams_number);
    std::random_device dev0;
    std::mt19937 rng0(dev0());
    double sigma_mirr=std::atof(argv[2]);   
    std::normal_distribution<double> dist0{0, sigma_mirr};
    // std::string moves_file = "vert_moves_for_"+std::to_string(2*beams_number)+"_beams.txt";
    std::string moves_file = "horiz_moves_for_"+std::to_string(2*beams_number)+"_beams.txt";
    // std::string moves_file = "angz_moves_for_"+std::to_string(2*beams_number)+"_beams.txt";
    // std::string moves_file = "angy_moves_for_"+std::to_string(2*beams_number)+"_beams.txt";
    bool rread=false;

    std::ifstream inm(moves_file);
    std::string line, firstline;
    size_t ii=0;   
    if (inm.is_open()){
        rread=true;
        getline (inm,firstline);
        while (getline (inm,line)) {
            mirr_adj[ii]=std::stod(line);
            ii++;
        }
        inm.close();
    } else {
        std::cout << "Unable to open file" << '\n';
    }


    std::ofstream outm(moves_file);
    if (rread)
    {
        std::cout << "meow1" << std::endl;
        outm << firstline << '\n';
        for (size_t i=0;i<2*beams_number;i++){
            outm << mirr_adj[i] << '\n';
        } 
    } else {
        std::cout << "meow2" << std::endl;
        outm << "for sigma = " <<  sigma_mirr << '\n';
        for (size_t i=0;i<2*beams_number;i++){
            mirr_adj[i]=dist0(rng0);
            outm << mirr_adj[i] << '\n';
        }
    }
    outm.close();

    double theta0;

    for (size_t i=0;i<beams_number;i++){


        // //double belt, angles_z
        // axFm1={F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)};
        // ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i), 1.0*cos(gamma_med)};
        // ax2m1=(ax1m1 % axFm1);
        // ax2m1 /= ax2m1.norm();
        // ParabolicSurface mirrorold1(
        // {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)},
        // ax1m1,ax2m1, 4.0*F, 4.0*F
        // );   

        // theta0=gamma_med+mirr_adj[i];
        // axFm1={F*cos(theta0)*cos(2*M_PI/beams_number*i), F*cos(theta0)*sin(2*M_PI/beams_number*i), -F*sin(theta0)};
        // ax1m1={1.0*sin(theta0)*cos(2*M_PI/beams_number*i), 1.0*sin(theta0)*sin(2*M_PI/beams_number*i), 1.0*cos(theta0)};
        // ax2m1=(ax1m1 % axFm1);
        // ax2m1 /= ax2m1.norm();
        // ParabolicSurface mirrortilted1(
        // {F*cos(theta0)*cos(2*M_PI/beams_number*i), F*cos(theta0)*sin(2*M_PI/beams_number*i), -F*sin(theta0)},
        // ax1m1,ax2m1, 4.0*F, 4.0*F
        // );

        // Position tilt_r01={F*cos(theta0)*cos(2*M_PI/beams_number*i), F*cos(theta0)*sin(2*M_PI/beams_number*i), -F*sin(theta0)};
        // Position new_r01=tilt_r01+mirrorold1.point({p,0})-mirrortilted1.point({p,0});

        // mirrors1.emplace_back(ParabolicSurface(
        //     new_r01,
        //     ax1m1,ax2m1, 4.0*F, 4.0*F
        // ));


        // axFm2={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
        // ax1m2={-1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), -1.0*cos(gamma_med)};
        // ax2m2=(ax1m2 % axFm2);
        // ax2m2 /= ax2m2.norm();
        // ParabolicSurface mirrorold2(
        // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
        // ax1m2, ax2m2, 4.0*F, 4.0*F
        // );

        // theta0=gamma_med+mirr_adj[i+beams_number];
        // axFm2={-F*cos(theta0)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(theta0)*sin(2*M_PI/beams_number*i+add_angle), F*sin(theta0)};
        // ax1m2={-1.0*sin(theta0)*cos(2*M_PI/beams_number*i+add_angle), -1.0*sin(theta0)*sin(2*M_PI/beams_number*i+add_angle), -1.0*cos(theta0)};
        // ax2m2=(ax1m2 % axFm2);
        // ax2m2 /= ax2m2.norm();
        // ParabolicSurface mirrortilted2(
        // {-F*cos(theta0)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(theta0)*sin(2*M_PI/beams_number*i+add_angle), F*sin(theta0)},
        // ax1m2, ax2m2, 4.0*F, 4.0*F
        // );

        // Position tilt_r02={-F*cos(theta0)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(theta0)*sin(2*M_PI/beams_number*i+add_angle), F*sin(theta0)};
        // Position new_r02=tilt_r02+mirrorold2.point({p,0})-mirrortilted2.point({p,0});

        // mirrors2.emplace_back(ParabolicSurface(
        //     new_r02,
        //     ax1m2,ax2m2, 4.0*F, 4.0*F
        // ));
        
        // axFm1={F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)};
        // ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i), 1.0*cos(gamma_med)};
        // ax2m1=(ax1m1 % axFm1);
        // ax2m1 /= ax2m1.norm();
        // axFm2={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
        // ax1m2={-1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), -1.0*cos(gamma_med)};
        // ax2m2=(ax1m2 % axFm2);
        // ax2m2 /= ax2m2.norm();




        // //double belt, angles_y
        // double ax1m1_x=1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i);
        // double ax1m1_y=1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i);
        // double ax1m1_z=1.0*cos(gamma_med);
        // axFm1={F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)};
        // ax1m1={ax1m1_x, ax1m1_y, ax1m1_z};
        // ax2m1=(ax1m1 % axFm1);
        // ax2m1 /= ax2m1.norm();
        // ParabolicSurface mirrorold1(
        // {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)},
        // ax1m1,ax2m1, 4.0*F, 4.0*F
        // );   

        // // Position fornewax1=mirrorold1.point({p+beam_radius,0})-mirrorold1.point({p,0});
        // // Vector newax1={fornewax1[0],fornewax1[1],fornewax1[2]};
        // // newax1 /= newax1.norm();
        // theta0=mirr_adj[i];
        // Vector axFm1_new;
        // double M1[3][3] = {
        //     {cos(theta0)+(1-cos(theta0))*ax1m1_x*ax1m1_x,(1-cos(theta0))*ax1m1_x*ax1m1_y-sin(theta0)*ax1m1_z,(1-cos(theta0))*ax1m1_x*ax1m1_z+sin(theta0)*ax1m1_y},
        //     {(1-cos(theta0))*ax1m1_x*ax1m1_y+sin(theta0)*ax1m1_z, cos(theta0)+(1-cos(theta0))*ax1m1_y*ax1m1_y,(1-cos(theta0))*ax1m1_y*ax1m1_z-sin(theta0)*ax1m1_x},
        //     {(1-cos(theta0))*ax1m1_x*ax1m1_z-sin(theta0)*ax1m1_y,(1-cos(theta0))*ax1m1_y*ax1m1_z+sin(theta0)*ax1m1_x,cos(theta0)+(1-cos(theta0))*ax1m1_z*ax1m1_z}
        // };

        // for (size_t j=0; j<3;j++){
        //     axFm1_new[j]=M1[j][0]*axFm1[0]+M1[j][1]*axFm1[1]+M1[j][2]*axFm1[2];
        //     // ax1m1_new[j]=M1[j][0]*ax1m1[0]+M1[j][1]*ax1m1[1]+M1[j][2]*ax1m1[2];
        // }
        // Vector ax2m1_new;
        // ax2m1_new=(ax1m1 % axFm1_new);
        // ax2m1_new /= ax2m1_new.norm();
        // ParabolicSurface mirrortilted1(
        //     {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)},
        //     ax1m1,ax2m1_new, 4.0*F, 4.0*F
        // );

        // Position old_r01={F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)};
        // Position new_r01=old_r01+mirrorold1.point({p,0})-mirrortilted1.point({p,0});

        // mirrors1.emplace_back(ParabolicSurface(
        //     new_r01,
        //     ax1m1,ax2m1_new, 4.0*F, 4.0*F
        // ));

        // // mirrors1.emplace_back(ParabolicSurface(
        // //     {F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)},
        // //     ax1m1,ax2m1_new, 4.0*F, 4.0*F
        // // ));


        // double ax1m2_x=-1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle);
        // double ax1m2_y=-1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle);
        // double ax1m2_z=-1.0*cos(gamma_med);
        // axFm2={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
        // ax1m2={ax1m2_x, ax1m2_y, ax1m2_z};
        // ax2m2=(ax1m2 % axFm2);
        // ax2m2 /= ax2m2.norm();
        // ParabolicSurface mirrorold2(
        // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
        // ax1m2, ax2m2, 4.0*F, 4.0*F
        // );

        // // Position fornewax2=mirrorold2.point({p+beam_radius,0})-mirrorold2.point({p,0});
        // // Vector newax2={fornewax2[0],fornewax2[1],fornewax2[2]};
        // // newax2 /= newax2.norm();
        // theta0=mirr_adj[i+beams_number];
        // Vector axFm2_new;
        // double M2[3][3] = {
        //     {cos(theta0)+(1-cos(theta0))*ax1m2_x*ax1m2_x,(1-cos(theta0))*ax1m2_x*ax1m2_y-sin(theta0)*ax1m2_z,(1-cos(theta0))*ax1m2_x*ax1m2_z+sin(theta0)*ax1m2_y},
        //     {(1-cos(theta0))*ax1m2_x*ax1m2_y+sin(theta0)*ax1m2_z, cos(theta0)+(1-cos(theta0))*ax1m2_y*ax1m2_y,(1-cos(theta0))*ax1m2_y*ax1m2_z-sin(theta0)*ax1m2_x},
        //     {(1-cos(theta0))*ax1m2_x*ax1m2_z-sin(theta0)*ax1m2_y,(1-cos(theta0))*ax1m2_y*ax1m2_z+sin(theta0)*ax1m2_x,cos(theta0)+(1-cos(theta0))*ax1m2_z*ax1m2_z}
        // };

        // for (size_t j=0; j<3;j++){
        //     axFm2_new[j]=M2[j][0]*axFm2[0]+M2[j][1]*axFm2[1]+M2[j][2]*axFm2[2];
        //     // ax1m2_new[j]=M1[j][0]*ax1m2[0]+M1[j][1]*ax1m2[1]+M1[j][2]*ax1m2[2];
        // }
        // Vector ax2m2_new;
        // ax2m2_new=(ax1m2 % axFm2_new);
        // ax2m2_new /= ax2m2_new.norm();
        // ParabolicSurface mirrortilted2(
        // {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
        // ax1m2, ax2m2_new, 4.0*F, 4.0*F
        // );

        // Position old_r02={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
        // Position new_r02=old_r02+mirrorold2.point({p,0})-mirrortilted2.point({p,0});

        // mirrors2.emplace_back(ParabolicSurface(
        //     new_r02,
        //     ax1m2, ax2m2_new, 4.0*F, 4.0*F
        // ));

        // // mirrors2.emplace_back(ParabolicSurface(
        // //     {-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)},
        // //     ax1m2, ax2m2_new, 4.0*F, 4.0*F
        // // ));



        //double belt, vert+horiz
        axFm1={F*cos(gamma_med)*cos(2*M_PI/beams_number*i), F*cos(gamma_med)*sin(2*M_PI/beams_number*i), -F*sin(gamma_med)};
        ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i), 1.0*cos(gamma_med)};
        ax2m1=(ax1m1 % axFm1);
        ax2m1 /= ax2m1.norm();
        mirrors1.emplace_back(ParabolicSurface(
            {(F+mirr_adj[i]/cos(gamma_med))*cos(gamma_med)*cos(2*M_PI/beams_number*i), (F+mirr_adj[i]/cos(gamma_med))*cos(gamma_med)*sin(2*M_PI/beams_number*i), -(F/*+mirr_adj[i]*/)*sin(gamma_med)/*+mirr_adj[i]*/},
            ax1m1,ax2m1, 4.0*F, 4.0*F
        ));

        axFm2={-F*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -F*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), F*sin(gamma_med)};
        ax1m2={-1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), -1.0*cos(gamma_med)};
        ax2m2=(ax1m2 % axFm2);
        ax2m2 /= ax2m2.norm();
        mirrors2.emplace_back(ParabolicSurface(
            {-(F+mirr_adj[i+beams_number]/cos(gamma_med))*cos(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), -(F+mirr_adj[i+beams_number]/cos(gamma_med))*cos(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), (F/*+mirr_adj[i+beams_number]*/)*sin(gamma_med)/*-mirr_adj[i+beams_number]*/},
            ax1m2, ax2m2, 4.0*F, 4.0*F
        ));


        beams1.emplace_back(ParallelBeamAlpha(lambda, {0.0, 0.0, 0.0}, ax1m1, -ax2m1,
                                //  [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 - h) + sqr(x2)) / (2 * sqr(sigma)) ); },
                                [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 - h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
                                //[h](double x1, double x2) { return sqr(x1 - h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
                                // [h](double x1, double x2) { return smoothed(sqrt(sqr(x1 - h) + sqr(x2)), beam_radius); },
                                [](double, double) { return 0.0; }));

        ax1m1={1.0*sin(gamma_med)*cos(2*M_PI/beams_number*i+add_angle), 1.0*sin(gamma_med)*sin(2*M_PI/beams_number*i+add_angle), 1.0*cos(gamma_med)};

        beams2.emplace_back(ParallelBeamAlpha(lambda, {0.0, 0.0, 0.0}, ax1m1, ax2m2,
                                //  [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 + h) + sqr(x2)) / (2 * sqr(sigma)) ); },
                                [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (pow(x1 + h,12) + pow(x2,12)) / (pow(sigma, 12)) ); },
                                //[h](double x1, double x2) { return sqr(x1 + h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
                                // [h](double x1, double x2) { return smoothed(sqrt(sqr(x1 + h) + sqr(x2)), beam_radius); },
                                [](double, double) { return 0.0; }));

        // refs1.emplace_back(mirrors1[i], beams1[i], region1);
        // refs2.emplace_back(mirrors2[i], beams2[i], region2);
        std::string filename_suffix = std::to_string(i);

        //plotting mirrors
        // plot_field_on_given_surface(mirrors1[i], beams1[i], region1, 100, "E", "mirror1", filename_suffix);
        // std::cout << "mirror1 plotted" << std::endl;
        // plot_field_on_given_surface(mirrors2[i], beams2[i], region2, 100, "E", "mirror2", filename_suffix);
        // std::cout << "mirror2 plotted" << std::endl;
    
    }
    std::cout << "mirror plotting ended" << std::endl;
    std::vector<StrattonChuReflection> refs;
    //double belt
    refs.reserve(2*beams_number);
    for (size_t i=0;i<2*beams_number;i++){
        if (i<beams_number){
            refs.emplace_back(mirrors1[i], beams1[i], region1);
        }
        else {
            refs.emplace_back(mirrors2[i-beams_number], beams2[i-beams_number], region2);
        }        
    }

    // //single belt
    // refs.reserve(beams_number);
    // for (size_t i=0;i<beams_number;i++){
    //     refs.emplace_back(mirrors1[i], beams1[i], region1);              
    // }
 


    // PlaneSurface focal_plane( p_focus, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} );
    // SurfaceRegion focal_region;
    // // focal_region.x_min = 0.0 - 12 * lambda;
    // // focal_region.x_max = 0.0 + 12 * lambda;
    // // focal_region.y_min = 0.0 - 12 * lambda;
    // // focal_region.y_max = 0.0 + 12 * lambda;    
    // focal_region.x_min = 0.0 - 3 * lambda;
    // focal_region.x_max = 0.0 + 3 * lambda;
    // focal_region.y_min = 0.0 - 3 * lambda;
    // focal_region.y_max = 0.0 + 3 * lambda;

    // PlaneSurface focal_plane_transversal( p_focus, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0} );
    // SurfaceRegion focal_region_transversal;
    // focal_region_transversal.x_min = 0.0 - 10 * lambda;
    // focal_region_transversal.x_max = 0.0 + 10 * lambda;
    // focal_region_transversal.y_min = 0.0 - 10 * lambda;
    // focal_region_transversal.y_max = 0.0 + 10 * lambda;

    PrlppdVolume focal_volume( p_focus, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
    VolumeRegion focal_area;
    double xcent=0.0;
    // double xcent=8 * lambda;
    // double xcent=-8 * lambda;
    double ycent=0.0;
    // double ycent=8 * lambda;
    // double ycent=-8 * lambda;
    double zcent=0.0;
    // double zcent=8 * lambda;
    // double zcent=-8 * lambda;
    focal_area.x_min = xcent - 4 * lambda;
    focal_area.x_max = xcent + 4 * lambda;
    focal_area.y_min = ycent - 4 * lambda;
    focal_area.y_max = ycent + 4 * lambda;
    focal_area.z_min = zcent - 4 * lambda;
    focal_area.z_max = zcent + 4 * lambda;

    std::cout << "meow" << std::endl;
    //double belt


    std::string nameE=std::string(GET_VARIABLE_NAME(focal_area))+"_"+std::to_string(2*beams_number)+"-beams_new_config";
    // std::cout << nameE << std::endl;
    FieldsIn3dRegion <focal_points> fields(focal_volume, focal_area, nameE, refs, 2*beams_number);    
    std::cout <<"all area ready" << std::endl;

    // std::string name=std::string(GET_VARIABLE_NAME(focal_region))+"_"+std::to_string(2*beams_number)+"-beams_new_config2";
    // FieldsInRegion <focal_points> fields_long(focal_plane, focal_region, name, refs, 2*beams_number);    
    // std::cout <<"fp ready" << std::endl;
    // name=std::string(GET_VARIABLE_NAME(focal_region_transversal))+"_"+std::to_string(2*beams_number)+"-beams_new_config2";
    // FieldsInRegion <focal_points_transversal> fields_trans(focal_plane_transversal, focal_region_transversal, name, refs, 2*beams_number);
    // std::cout <<"fp trans ready" << std::endl;

    // //single belt
    // std::string name=std::string(GET_VARIABLE_NAME(focal_region))+"_"+std::to_string(beams_number)+"-beams_single";
    // FieldsInRegion <focal_points> fields_long(focal_plane, focal_region, name, refs, beams_number);    
    // std::cout <<"fp ready" << std::endl;
    // name=std::string(GET_VARIABLE_NAME(focal_region_transversal))+"_"+std::to_string(beams_number)+"-beams_single";
    // FieldsInRegion <focal_points_transversal> fields_trans(focal_plane_transversal, focal_region_transversal, name, refs, beams_number);
    // std::cout <<"fp trans ready" << std::endl;

    // std::vector<double> phase_lims(6);
    // phase_lims[0]=M_PI/3;
    // phase_lims[1]=M_PI/4;
    // phase_lims[2]=M_PI/6;
    // phase_lims[3]=0;
    // phase_lims[4]=M_PI/2;
    // phase_lims[5]=M_PI;
    std::vector<double> phase_lims(5);
    phase_lims[0]=0;
    phase_lims[1]=M_PI/3;
    phase_lims[2]=M_PI/2;
    phase_lims[3]=M_PI;
    phase_lims[4]=2*M_PI;

    std::vector<double> phases;
    phases.reserve(2*beams_number);
    // phases.reserve(beams_number);

    std::random_device dev;
    std::mt19937 rng(dev());

    double testIntens[5][10];
    double trueint;

    std::string intens_file = "intens_for_"+std::to_string(2*beams_number)+"_beams.txt";
    // std::string intens_file = "intens_for_"+std::to_string(beams_number)+"_beams_single.txt";
    std::ofstream outr(intens_file);


    // std::string filename_suffix1;
    // for (size_t i=0;i<2*beams_number;i++){
    //     if (i<beams_number){
    //         filename_suffix1 = "_beam"+std::to_string(i+1)+"_side1";
    //     }
    //     else {
    //         filename_suffix1 = "_beam"+std::to_string(i+1-beams_number)+"_side2";
    //     }    
    //     // plot_field_on_given_surface(planes[i], refs[i], focal_region, focal_points, "E1", "longitudinal_E1", filename_suffix1);
    //     plot_field_on_given_surface_3darr(planes[i], focal_region, fields_long, 0, focal_points, "E1", "longitudinal_E1", filename_suffix1);
    //     // plot_field_on_given_surface_3darr(focal_plane, focal_region, fields_long, i, focal_points, "E1", "longitudinal_E1", filename_suffix1);
    //     // plot_field_on_given_surface_3darr(focal_plane_transversal, focal_region_transversal, fields_trans, i, focal_points_transversal, "E1", "transversal_E1", filename_suffix1);
    //     std::cout << "meow" << std::endl;
    // }
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
                // std::cout << phases[i] << std::endl;
            }
            std::string filename_suffix = "_"+std::to_string(2*beams_number)+"beams_test"+std::to_string(phase_iter)+"_"+std::to_string(test_iter);
            // std::string filename_suffix = "_"+std::to_string(beams_number)+"beams_single_test"+std::to_string(phase_iter)+"_"+std::to_string(test_iter);

            // fields_long.constructSummaryMono(phases);
            // plot_field_on_given_surface_3darr(focal_plane, focal_region, fields_long, focal_points, "E1", "longitudinal_E1", filename_suffix);
            // fields_trans.constructSummaryMono(phases);
            // plot_field_on_given_surface_3darr(focal_plane_transversal, focal_region_transversal, fields_trans, focal_points_transversal, "E1", "transversal_E1", filename_suffix);
            fields.constructSummaryMono(phases);
            // plot_field_in_given_area_3darr(focal_volume, focal_area, fields, focal_points, "E1", "field", filename_suffix);

            // fields_long.findSummaryMax();
            // fields_long.findSummaryMaxExact();
            // if (fields_long.checkMax()){
            //     std::cout <<"ok"<< std::endl;
            // }
            // else{
            //     std::cout <<"not ok"<< std::endl;
            //     trueint=fields_long.calculateTrueMaxIntensitySummary();
            //     std::cout <<"test "<< test_iter <<" max intens exact = " << trueint << std::endl;
            // }
            // testIntens[phase_iter][test_iter]=fields_long.calculateMaxIntensitySummary();
            // std::cout <<"test "<< test_iter <<" max intens = " << testIntens[phase_iter][test_iter] << std::endl;
            // Vector2D MaxPos=fields_long.MaxIntensityPoint();            
            // outr << testIntens[phase_iter][test_iter] << "   in point [" << MaxPos[0] << ", " << MaxPos[1] <<"]" <<'\n';

            fields.findSummaryMax();
            fields.findSummaryMaxExact();
            if (fields.checkMax()){
                std::cout <<"ok"<< std::endl;
            }
            else{
                std::cout <<"not ok"<< std::endl;
                trueint=fields.calculateTrueMaxIntensitySummary();
                std::cout <<"test "<< test_iter <<" max intens exact = " << trueint << std::endl;
            }
            testIntens[phase_iter][test_iter]=fields.calculateMaxIntensitySummary();
            std::cout <<"test "<< test_iter <<" max intens = " << testIntens[phase_iter][test_iter] << std::endl;
            Vector MaxPos=fields.MaxIntensityPoint();
            
            outr << testIntens[phase_iter][test_iter] << "   in point [" << MaxPos[0] << ", " << MaxPos[1] << ", " << MaxPos[2] <<"]" <<'\n';

        // }
    }
    
    outr.close();


    //конец этой части кода

//     std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

//     ///немонохроматический пучок
//     const double phi = M_PI / 3;
//     double alpha = M_PI / 6;
//     const double omega0 = 2 * M_PI / lambda*c;   
    
//     double tauint; //для интенсивности
//     tauint=11e-15;
//     //std::cout << tauint << std::endl;

//     double tauf=tauint/(2*sqrt(2*log(2))); //для поля
//     //std::cout << tauf << std::endl;

//     //funct=exp(-sqr(t) / (4 * sqr(tauf)) * exp(-Complex(0.0, 1.0) * omega0 * t);

//     fftw_complex *in, *out;
//     fftw_plan plan;
//     double T=100*tauf;
//     int numb=5000;
//     //int numb=625;
//     //int numb=100;
//     double deltat=T/(numb-1);
//     double tstart=-T/10;
//     //reinterpret_cast //для совместимости форматов комплексных чисел
//     std::complex<double> funct;


//     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*numb);
//     for (int i = 0; i < numb; i++)
//     {
//         funct=exp(-sqr(tstart+deltat*i) / (4 * sqr(tauf))) * exp(std::complex<double>(0.0, 1.0) * omega0 * (tstart+deltat*i));
//         //funct=exp(-sqr(tstart+deltat*i) / (4 * sqr(tauf))) * cos(omega0 * (tstart+deltat*i));
//         in[i][0]=funct.real();
//         in[i][1]=funct.imag();
//     }
//     /*
//     ///тест
//     std::ofstream inr("input.txt");
// 	for (int i = 0; i < numb; i++) {
//         if (in[i][1]<0)
//             inr << in[i][0] <<"-"<<abs(in[i][1])<<"i"<< '\n';
//         else
//             inr << in[i][0] <<"+"<<in[i][1]<<"i"<< '\n';
// 	}
// 	inr.close();

//     std::ofstream iuto("inputabs.txt");
// 	for (int i = 0; i < numb; i++) {
//          iuto << sqrt(sqr(in[i][0])+sqr(in[i][1])) << '\n';
// 	}
// 	iuto.close();
//     */
//     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*numb);
//     plan = fftw_plan_dft_1d (numb, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

//     fftw_execute (plan);

//     std::complex<double> *outc=reinterpret_cast<std::complex<double>*>(out);

//     // //тест
//     // std::ofstream outr("output.txt");
// 	// for (int i = 0; i < numb; i++) {
//     //     if (out[i][1]<0)
//     //         outr << out[i][0] <<"-"<<abs(out[i][1])<<"i"<< '\n';
//     //     else
//     //         outr << out[i][0] <<"+"<<out[i][1]<<"i"<< '\n';
// 	// }
// 	// outr.close();

//     // std::ofstream outre("output1.txt");
// 	// for (int i = 0; i < numb; i++) {
//     //     outre << outc[i] << '\n';
// 	// }
// 	// outre.close();

//     std::ofstream outo("outputabs.txt");
// 	for (int i = 0; i < numb; i++) {
//          //outo << sqrt(sqr(out[i][0])+sqr(out[i][1])) << '\n';
//          outo << abs(outc[i]) << '\n';
//          //std::cout << sqrt(sqr(out[i][0])+sqr(out[i][1])) << std::endl;
// 	}
// 	outo.close();



//     int central_harmonic=round(T*omega0/2/M_PI);
//     double halfshir=M_PI/tauf;
//     int left_harmonics_border=round(T*(omega0-halfshir)/2/M_PI);
//     int right_harmonics_border=round(T*(omega0+halfshir)/2/M_PI);

//     ///исправления для учёта всех гармоник
//     int harmonics_number_full=right_harmonics_border-left_harmonics_border;
//     std::cout << harmonics_number_full << std::endl;
//     if (harmonics_number_full==harmonics_number) {
//         std::cout << "Calculations are made for all points in spectrum" << std::endl;
//     } else {
//         std::cout << "Set harmonics_number in constants.hpp to " << harmonics_number_full << std::endl;
//     }

//     ///для учёта избранного числа гармоник
//     // int harmonic_interval=round(static_cast<double>(right_harmonics_border-left_harmonics_border)/(harmonics_number+1));
//     // int middle_index=harmonics_number/2;




//     //std::cout << central_harmonic << std::endl;


//        // Параметры для гауссова пучка, той же дисперсии, что и _|¯¯¯|_, и чтобы в нём энергии как в _|¯¯¯|_
//     double sigma = beam_radius / 2;
//     //double gauss_A = 2/sqrt(harmonics_number);
//     double gauss_A = 2;
//     //std::vector<double> gauss_AA(harmonics_number);

//     ///для учёта избранного числа гармоник
//     // std::vector<double> freqs(harmonics_number);
//     // std::vector<double> lambdas(harmonics_number);
//     // std::vector<std::complex<double>> ampls(harmonics_number);
//     // freqs[middle_index]=omega0;
//     // lambdas[middle_index]=lambda;
//     // ampls[middle_index]=outc[central_harmonic];

//     //gauss_AA[middle_index]=2/sqrt(harmonics_number);

//     ///исправления для учёта всех гармоник
//     std::vector<double> freqs(harmonics_number_full);
//     std::vector<double> lambdas(harmonics_number_full); 
//     std::vector<std::complex<double>> ampls(harmonics_number_full);
//     // std::vector<double> freqs_full(harmonics_number_full);
//     // std::vector<double> lambdas_full(harmonics_number_full); 
//     // std::vector<std::complex<double>> ampls_full(harmonics_number_full);

//     for (int i = 0; i < harmonics_number_full; i++) {
//         freqs[i]=(left_harmonics_border+i)/T*2*M_PI;
//         lambdas[i]=2*M_PI*c/freqs[i];
//         ampls[i]=outc[left_harmonics_border+i];
//         //std::cout << freqs[i] << std::endl;
//     }

//     freqs[central_harmonic-left_harmonics_border]=omega0;
//     lambdas[central_harmonic-left_harmonics_border]=lambda;

// /*
//     if (middle_index!=0) {
//         for (int i = 1; i <= middle_index; i++) {
//             freqs[middle_index+i]=(central_harmonic+harmonic_interval*i)/T*2*M_PI;
//             lambdas[middle_index+i]=2*M_PI*c/freqs[middle_index+i];
//             ampls[middle_index+i]=outc[central_harmonic+harmonic_interval*i];
//             freqs[middle_index-i]=(central_harmonic-harmonic_interval*i)/T*2*M_PI;
//             lambdas[middle_index-i]=2*M_PI*c/freqs[middle_index-i];
//             ampls[middle_index-i]=outc[central_harmonic-harmonic_interval*i];
            
//             //std::cout << freqs[i] << std::endl;

//             //gauss_AA[middle_index+i]=2/sqrt(harmonics_number)*abs(outc[central_harmonic])/abs(outc[central_harmonic+harmonic_interval*i]);
//             //gauss_AA[middle_index-i]=2/sqrt(harmonics_number)*abs(outc[central_harmonic])/abs(outc[central_harmonic-harmonic_interval*i]);
//         }
//     }
// */


//     // for (int i = 0; i < harmonics_number; i++) {
//     //         std::cout << gauss_AA[i] << std::endl;
//     // }
    
    

//     double F = get_F_by_beam_parameters_alpha(alpha, phi, beam_width);
//     double p = get_p_by_beam_parameters_alpha(alpha, F); // impact parameter по нижней границе
//     double h = p + beam_radius;

//     std::cout << "\t Focal length = " << F<< std::endl;
//     /*if (plot_distortion) {std::cout << "\t k = " << distortion_k << "\t ampl = " << distortion_ampl;}
//     std::cout << std::endl;*/

//     Position p_focus = {0.0, 0.0, 0.0};

//     ParabolicSurface mirror1(
//         {0.0, 0.0, -F},
//         {1.0, 0.0, 0.0},
//         {0.0, 1.0, 0.0},
//         4.0*F, 4.0*F
//     );

//     ParabolicSurface mirror2(
//         {0.0, 0.0, F},
//         {-1.0, 0.0, 0.0},
//         {0.0, 1.0, 0.0},
//         4.0*F, 4.0*F
//     );

//     PlaneSurface beam_profile({0.0, 0.0, 30.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0});


//     //std::cout << gauss_A << std::endl;

//     double mirror_radius = 3 * sigma;
//     //double mirror_radius = beam_radius;
//     //double mirror_radius = beam_radius + 0.1;

//     SurfaceRegion region1;
//     region1.x_min = h - mirror_radius;
//     region1.x_max = h + mirror_radius;
//     region1.y_min = 0.0 - mirror_radius;
//     region1.y_max = 0.0 + mirror_radius;

//     SurfaceRegion region2; //           Почему не минус h ???
//     region2.x_min = h - mirror_radius;
//     region2.x_max = h + mirror_radius;
//     region2.y_min = 0.0 - mirror_radius;
//     region2.y_max = 0.0 + mirror_radius;

//     SurfaceRegion region_profile;
//     region_profile.x_min = h - 1.1 * mirror_radius;
//     region_profile.x_max = h + 1.1 * mirror_radius;
//     region_profile.y_min = 0.0 - 1.1 * mirror_radius;
//     region_profile.y_max = 0.0 + 1.1 * mirror_radius;

//     std::vector<ParallelBeamAlpha> beams;
//     std::vector<StrattonChuReflection> refs;
//     beams.reserve(harmonics_number);
//     refs.reserve(harmonics_number);
//     for (size_t i=0;i<harmonics_number;i++)
//     {
//         //double gauss_A=gauss_AA[i];
//         //std::cout << gauss_A << "  "<< gauss_AA[i]<< std::endl;
//          beams.emplace_back(ParallelBeamAlpha(lambdas[i], {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0},
//                             [h, gauss_A, sigma](double x1, double x2) { return gauss_A * exp( - (sqr(x1 - h) + sqr(x2)) / (2 * sqr(sigma)) ); },
//                             //[h](double x1, double x2) { return sqr(x1 - h) + sqr(x2) <= sqr(beam_radius) ? 1.0 : 0.0; },
//                             //[h](double x1, double x2) { return smoothed(sqrt(sqr(x1 - h) + sqr(x2)), beam_radius); },
//                             [](double, double) { return 0.0; }));
//          refs.emplace_back(mirror1, beams[i], region1);
//     }

//     /*
//     std::cout << refs.size() << std::endl;
//     std::cout << (sizeof(freqs)/sizeof(*freqs)) << std::endl;
//     std::cout << (sizeof(ampls)/sizeof(*ampls)) << std::endl;
//     */

//     //cache method
//     //SpecIFT impulse(refs, freqs, ampls);

    
//     PlaneSurface focal_plane( p_focus, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0} );
//     SurfaceRegion focal_region;
//     focal_region.x_min = 0.0 - 12 * lambda;
//     focal_region.x_max = 0.0 + 12 * lambda;
//     focal_region.y_min = 0.0 - 12 * lambda;
//     focal_region.y_max = 0.0 + 12 * lambda;


//     Vector direction_vector = p_focus - mirror1.point({h, 0.0});
//     direction_vector /= direction_vector.norm();

//     PlaneSurface focal_plane_transversal( p_focus, {0.0, 1.0, 0.0}, Vector(0.0, 1.0, 0.0) % direction_vector );
//     SurfaceRegion focal_region_transversal;
//     focal_region_transversal.x_min = 0.0 - 10 * lambda;
//     focal_region_transversal.x_max = 0.0 + 10 * lambda;
//     focal_region_transversal.y_min = 0.0 - 10 * lambda;
//     focal_region_transversal.y_max = 0.0 + 10 * lambda;


//     ///3d-array method
//     // std::string name1=std::string(GET_VARIABLE_NAME(focal_region))+"_tauint-"+std::to_string(tauint*1e15);
//     // std::string name2=std::string(GET_VARIABLE_NAME(focal_region_transversal))+"_tauint-"+std::to_string(tauint*1e15);
//     std::string name1=std::string(GET_VARIABLE_NAME(focal_region))+"_alpha-"+std::to_string(int(alpha/M_PI*180));
//     std::string name2=std::string(GET_VARIABLE_NAME(focal_region_transversal))+"_alpha-"+std::to_string(int(alpha/M_PI*180));
    
//     // FieldsInRegion <focal_points> pulse1long(focal_plane, WITH_NAME(focal_region), refs, harmonics_number);    
//     FieldsInRegion <focal_points> pulse1long(focal_plane, focal_region, name1, refs, harmonics_number);    
//     std::cout <<"fp ready" << std::endl;
//     FieldsInRegion <focal_points_transversal> pulse1trans(focal_plane_transversal, focal_region_transversal, name2, refs, harmonics_number);
//     std::cout <<"fp trans ready" << std::endl;
    
    

//     double Tmax=5*tauf;
//     int number_of_time_points=25;
//     double time_interval=Tmax/number_of_time_points;


//     std::cout << harmonics_number << " harmonics modeling starts" << std::endl;

//     ///cache method
//     for (size_t i=0;i<number_of_time_points;i++)
//     {   
//         std::string filename_suffix = "time-";
//         filename_suffix += std::to_string(i);
//         std::cout <<"file" << i  << std::endl;
//         //std::cout << "freqs "<< (sizeof(freqs)/sizeof(*freqs)) << std::endl;
//         //std::cout << "ampls "<< (sizeof(ampls)/sizeof(*ampls)) << std::endl;

//         // //cache method
//         // plot_field_on_given_surface_with_time(focal_plane, impulse, focal_region, focal_points, i*time_interval, "E1", "longitudinal_E1", filename_suffix);
//         // std::cout << "longitudinal_E1 plotted" << std::endl;
    
//         // plot_field_on_given_surface_with_time(focal_plane_transversal, impulse, focal_region_transversal, focal_points_transversal, i*time_interval, "E1", "transversal_E1", filename_suffix);
//         // std::cout << "transversal_E1 plotted" << std::endl;
//         // if (i==0)  {
//         //     using CACHE=Caching::ConcurrentCache<Caching::Dependances<double,double,double,double>, FieldValue>;
//         //     auto& cache = STATIC_LOCAL_GET_REF(CACHE, conc_cache);
//         //     cache.set_stores_availability(false);
//         // }
//         //3d-array method
//         plot_field_on_given_surface_with_time_3darr(focal_plane, focal_region, pulse1long, freqs, ampls, focal_points, i*time_interval, "E1", "longitudinal_E1", filename_suffix);
//         std::cout << "longitudinal_E1 plotted" << std::endl;
    
//         plot_field_on_given_surface_with_time_3darr(focal_plane_transversal, focal_region_transversal, pulse1trans, freqs, ampls, focal_points_transversal, i*time_interval, "E1", "transversal_E1", filename_suffix);
//         std::cout << "transversal_E1 plotted" << std::endl;

//     }

    
//     std::cout << harmonics_number << " harmonics modeling is over" << std::endl;

//     fftw_destroy_plan(plan);
//     fftw_free(in);
//     fftw_free(out);

//     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//     std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000.0/60.0 << " min" << std::endl;
  
    return 0;

}
