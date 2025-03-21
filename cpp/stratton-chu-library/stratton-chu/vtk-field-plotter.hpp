#ifndef VTKFIELDPLOTTER_HPP
#define VTKFIELDPLOTTER_HPP

#include "stratton-chu/surface.hpp"
#include "stratton-chu/volume.hpp"

#include "stratton-chu/parallel-beam.hpp"
#include "stratton-chu/stratton-chu-field.hpp"
#include "stratton-chu/field.hpp"
#include "stratton-chu/utils.hpp"

#include "stratton-chu/csv-saver.hpp"
#include "stratton-chu/vtk-saver.hpp"

#include <iostream>
#include "threading.hpp"

#include <cmath>

#include <fstream>
#include <string>

#include <vector>

#include "cache.hpp"
#include "static_local_tracker.hpp"


void plot_field_on_given_surface(
        /*const*/ ISurface& surface,
        /*const*/ IField& field,
        /*const*/ SurfaceRegion& region,
        // double phase,
        int n_points,
        const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix
        )
{
    VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
    VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());
    VTKSurfaceSaver vtk_saver_int(n_points, n_points, (quantity_name + "_intens").c_str());

    // В тау-тау-эн;
    VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

    struct Task{
            SurfaceRegion &region;
            IField &field;
            ISurface &surface;
            VTKSurfaceSaver &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_int, &vtk_saver_surf_r, &vtk_saver_surf_m;
            int n_points;
            size_t  i, j;
            void operator()(){
                Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                Position p = surface.point(xy);
                // if (i==0) {std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;}
                // // std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;

                FieldValue field_value = field.get(p);
                // if (n_points==10) {std::cout << "point " << i << " " << j << " ready" << std::endl;}
                // if (n_points==10) {if (j==0) {std::cout << "point ready for 10" << std::endl;}}
                //std::cout << "point ready" << std::endl;
                // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

                vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
                vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
                vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

                double absE;
                Vector ves;
                absE=0.0;
                ves=max_field(field_value.E);
                for (size_t l=0;l<3;l++){
                    absE+=ves[l]*ves[l];
                }
                absE=absE*4e13*50/16;
                ves[0]=absE;
                ves[1]=0;
                ves[2]=0;
                vtk_saver_int.set_point(i, j, p, ves);

                FieldValue field_rotated;
                field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
                field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

                // field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
                // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                // field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

                // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
                // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));

                vtk_saver_surf_r.set_point(i, j, p, vec_modulus(field_value.E));
                vtk_saver_surf_m.set_point(i, j, p, vec_phases(field_value.E));
            }
            operator bool()const{
                return true;
            }
        };
    using Threads = Threading::ThreadPool<ConcurrentQueue<Task>,14>;
    typename Threads::TaskQueue tasks;
    // using Threads = ThreadPool<std::function<void()>>;
    // Threads::TaskQueue tasks;

    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            tasks.push(Task{region, field, surface, vtk_saver_r, vtk_saver_i, vtk_saver_m, vtk_saver_int, vtk_saver_surf_r, vtk_saver_surf_m, n_points, i, j});
            // {
            //     Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
            //     Position p = surface.point(xy);
            //     FieldValue field_value = field.get(p);
            //     //std::cout << "point ready" << std::endl;
            //     // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

            //     vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
            //     vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
            //     vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

            //     FieldValue field_rotated;
            //     // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
            //     // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
            //     // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

            //     field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
            //     field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
            //     field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

            //     vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
            //     vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
            // });
        }
    }
    Threads threads{tasks};
    threads.run();
    // Threads threads{tasks};
    // threads.run();

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());
    vtk_saver_int.save((filename_prefix + "_intens_" + filename_suffix).c_str());

    // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
    vtk_saver_surf_r.save((filename_prefix + "_modulus_" + filename_suffix).c_str());
    vtk_saver_surf_m.save((filename_prefix + "_phases_" + filename_suffix).c_str());
}


void plot_field_on_given_surface_ell(
        /*const*/ ISurface& surface,
        /*const*/ IField& field,
        /*const*/ SurfaceRegion& region,
        // double phase,
        int n_points,
        const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix
        )
{
    VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
    VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());
    VTKSurfaceSaver vtk_saver_int(n_points, n_points, (quantity_name + "_intens").c_str());

    // В тау-тау-эн;
    VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

    struct Task{
            SurfaceRegion &region;
            IField &field;
            ISurface &surface;
            VTKSurfaceSaver &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_int, &vtk_saver_surf_r, &vtk_saver_surf_m;
            int n_points;
            size_t  i, j;
            void operator()(){
                Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                Position p = surface.point(xy);
                // if (i==0) {std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;}
                // // std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;

                FieldValue field_value = field.get(p);
                std::cout << "point " << i << " " << j << " ready" << std::endl;
                // if (n_points==10) {if (j==0) {std::cout << "point ready for 10" << std::endl;}}
                //std::cout << "point ready" << std::endl;
                // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

                vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
                vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
                vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

                double absE;
                Vector ves;
                absE=0.0;
                ves=max_field(field_value.E);
                for (size_t l=0;l<3;l++){
                    absE+=ves[l]*ves[l];
                }
                absE=absE*4e13*50/16;
                ves[0]=absE;
                ves[1]=0;
                ves[2]=0;
                vtk_saver_int.set_point(i, j, p, ves);

                FieldValue field_rotated;
                field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
                field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

                // field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
                // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                // field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

                // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
                // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));

                vtk_saver_surf_r.set_point(i, j, p, vec_modulus(field_value.E));
                vtk_saver_surf_m.set_point(i, j, p, vec_phases(field_value.E));
            }
            operator bool()const{
                return true;
            }
        };
    using Threads = Threading::ThreadPool<ConcurrentQueue<Task>,14>;
    typename Threads::TaskQueue tasks;
    // using Threads = ThreadPool<std::function<void()>>;
    // Threads::TaskQueue tasks;

    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            tasks.push(Task{region, field, surface, vtk_saver_r, vtk_saver_i, vtk_saver_m, vtk_saver_int, vtk_saver_surf_r, vtk_saver_surf_m, n_points, i, j});
            // {
            //     Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
            //     Position p = surface.point(xy);
            //     FieldValue field_value = field.get(p);
            //     //std::cout << "point ready" << std::endl;
            //     // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

            //     vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
            //     vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
            //     vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

            //     FieldValue field_rotated;
            //     // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
            //     // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
            //     // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

            //     field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
            //     field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
            //     field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

            //     vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
            //     vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
            // });
        }
    }
    Threads threads{tasks};
    threads.run();
    // Threads threads{tasks};
    // threads.run();

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());
    vtk_saver_int.save((filename_prefix + "_intens_" + filename_suffix).c_str());

    // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
    vtk_saver_surf_r.save((filename_prefix + "_modulus_" + filename_suffix).c_str());
    vtk_saver_surf_m.save((filename_prefix + "_phases_" + filename_suffix).c_str());
}


void plot_field_on_given_surface_non_threaded(
        /*const*/ ISurface& surface,
        /*const*/ IField& field,
        /*const*/ SurfaceRegion& region,
        // double phase,
        int n_points,
        const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix
        )
{
    VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
    VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());

    // В тау-тау-эн;
    VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

    // struct Task{
    //         SurfaceRegion &region;
    //         IField &field;
    //         ISurface &surface;
    //         VTKSurfaceSaver &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m;
    //         int n_points;
    //         size_t  i, j;
    //         void operator()(){
    //             Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
    //             Position p = surface.point(xy);
    //             // if (i==0) {std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;}
    //             // // std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;

    //             FieldValue field_value = field.get(p);
    //             if (n_points==5) {if (j==0) {std::cout << "point ready" << std::endl;}}
    //             //std::cout << "point ready" << std::endl;
    //             // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

    //             vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
    //             vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
    //             vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

    //             FieldValue field_rotated;
    //             field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
    //             field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
    //             field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

    //             // field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
    //             // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
    //             // field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

    //             // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
    //             // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));

    //             vtk_saver_surf_r.set_point(i, j, p, vec_modulus(field_value.E));
    //             vtk_saver_surf_m.set_point(i, j, p, vec_phases(field_value.E));
    //         }
    //         operator bool()const{
    //             return true;
    //         }
    //     };
    // using Threads = Threading::ThreadPool<ConcurrentQueue<Task>,14>;
    // typename Threads::TaskQueue tasks;
    // using Threads = ThreadPool<std::function<void()>>;
    // Threads::TaskQueue tasks;

    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            // tasks.push(Task{region, field, surface, vtk_saver_r, vtk_saver_i, vtk_saver_m, vtk_saver_surf_r, vtk_saver_surf_m, n_points, i, j});
            Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
            Position p = surface.point(xy);
            // if (i==0) {std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;}
            // // std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;

            FieldValue field_value = field.get(p);
            if (n_points==10) {if (j==0) {std::cout << "point ready" << std::endl;}}
            //std::cout << "point ready" << std::endl;
            // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

            vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
            vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
            vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

            FieldValue field_rotated;
            field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
            field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
            field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

            // field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
            // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
            // field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

            // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
            // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));

            vtk_saver_surf_r.set_point(i, j, p, vec_modulus(field_value.E));
            vtk_saver_surf_m.set_point(i, j, p, vec_phases(field_value.E));
        }
    }
    // Threads threads{tasks};
    // threads.run();

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

    // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
    vtk_saver_surf_r.save((filename_prefix + "_modulus_" + filename_suffix).c_str());
    vtk_saver_surf_m.save((filename_prefix + "_phases_" + filename_suffix).c_str());
}

//smart - заполнитель массива с полями для последующего использования в интерполяции

template <const size_t array_points>
void plot_field_on_given_surface_smart(
        /*const*/ ISurface& surface,
        /*const*/ IField& field,
        /*const*/ SurfaceRegion& region,
        // double phase,
        int n_points,
        const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix,
        std::array<double, array_points> &xs, std::array<double, array_points> &ys, std::array<double, array_points> &zs, std::array<FieldValue, array_points> &FieldInRegion
        )
{
    VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
    VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());

    // В тау-тау-эн;
    VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

    struct Task{
            SurfaceRegion &region;
            IField &field;
            ISurface &surface;
            VTKSurfaceSaver &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m;
            int n_points;
            std::array<double, array_points> &xs;
            std::array<double, array_points> &ys;
            std::array<double, array_points> &zs;
            std::array<FieldValue, array_points*array_points> &FieldInRegion;
            size_t  i, j;
            void operator()(){
                Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                Position p = surface.point(xy);

                size_t count = i*n_points + j;

                xs[count]=p[0];
                ys[count]=p[1];
                zs[count]=p[2];

                // if (i==0) {std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;}
                // // std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;

                FieldValue field_value = field.get(p);
                FieldInRegion[count]=field_value;
                //std::cout << "point ready" << std::endl;
                // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);
                vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
                vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
                vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

                FieldValue field_rotated;
                field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
                field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

                // field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
                // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                // field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

                // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
                // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
                vtk_saver_surf_r.set_point(i, j, p, vec_modulus(field_value.E));
                vtk_saver_surf_m.set_point(i, j, p, vec_phases(field_value.E));
                
            }
            operator bool()const{
                return true;
            }
        };
    using Threads = Threading::ThreadPool<ConcurrentQueue<Task>,14>;
    typename Threads::TaskQueue tasks;
    // using Threads = ThreadPool<std::function<void()>>;
    // Threads::TaskQueue tasks;

    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            tasks.push(Task{region, field, surface, vtk_saver_r, vtk_saver_i, vtk_saver_m, vtk_saver_surf_r, vtk_saver_surf_m, n_points, xs, ys, zs, FieldInRegion, i, j});
            // {
            //     Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
            //     Position p = surface.point(xy);
            //     FieldValue field_value = field.get(p);
            //     //std::cout << "point ready" << std::endl;
            //     // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

            //     vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
            //     vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
            //     vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

            //     FieldValue field_rotated;
            //     // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
            //     // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
            //     // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

            //     field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
            //     field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
            //     field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

            //     vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
            //     vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
            // });
        }
    }
    Threads threads{tasks};
    threads.run();
    // Threads threads{tasks};
    // threads.run();

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

    // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
    vtk_saver_surf_r.save((filename_prefix + "_modulus_" + filename_suffix).c_str());
    vtk_saver_surf_m.save((filename_prefix + "_phases_" + filename_suffix).c_str());
}


template <const size_t array_points>
void plot_field_on_given_surface_smart_2d(
        /*const*/ ISurface& surface,
        /*const*/ IField& field,
        /*const*/ SurfaceRegion& region,
        // double phase,
        int n_points,
        const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix,
        std::array<double, array_points> &xs, std::array<double, array_points> &ys, std::array<std::array<FieldValue, array_points>, array_points> &FieldInRegion
        )
{
    VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
    VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());

    // В тау-тау-эн;
    VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());



    struct Task{
            SurfaceRegion &region;
            IField &field;
            ISurface &surface;
            VTKSurfaceSaver &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m;
            int n_points;
            std::array<double, array_points> &xs;
            std::array<double, array_points> &ys;
            std::array<std::array<FieldValue, array_points>, array_points> &FieldInRegion;
            size_t  i, j;
            void operator()(){
                Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                Position p = surface.point(xy);

                // size_t count = i*n_points + j;

                if (j==0) {xs[i]=xy[0];}
                if (i==0) {ys[j]=xy[1];}

                // if (i==0) {std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;}
                // // std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;

                FieldValue field_value = field.get(p);
                // FieldInRegion[count]=field_value;

                FieldInRegion[i][j]=field_value;

                //std::cout << "point ready" << std::endl;
                // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

                vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
                vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
                vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

                FieldValue field_rotated;
                // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
                // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

                field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
                field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

                // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
                // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));

                vtk_saver_surf_r.set_point(i, j, p, vec_modulus(field_value.E));
                vtk_saver_surf_m.set_point(i, j, p, vec_phases(field_value.E));
            }
            operator bool()const{
                return true;
            }
        };
    using Threads = Threading::ThreadPool<ConcurrentQueue<Task>,14>;
    typename Threads::TaskQueue tasks;
    // using Threads = ThreadPool<std::function<void()>>;
    // Threads::TaskQueue tasks;

    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            tasks.push(Task{region, field, surface, vtk_saver_r, vtk_saver_i, vtk_saver_m, vtk_saver_surf_r, vtk_saver_surf_m, n_points, xs, ys, FieldInRegion, i, j});
            // {
            //     Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
            //     Position p = surface.point(xy);
            //     FieldValue field_value = field.get(p);
            //     //std::cout << "point ready" << std::endl;
            //     // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

            //     vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
            //     vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
            //     vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

            //     FieldValue field_rotated;
            //     // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
            //     // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
            //     // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

            //     field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
            //     field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
            //     field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

            //     vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
            //     vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
            // });
        }
    }
    Threads threads{tasks};
    threads.run();
    // Threads threads{tasks};
    // threads.run();

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

    // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());

    vtk_saver_surf_r.save((filename_prefix + "_modulus_" + filename_suffix).c_str());
    vtk_saver_surf_m.save((filename_prefix + "_phases_" + filename_suffix).c_str());
}



// template <const size_t array_points>
// void plot_field_on_given_surface_smart_2d_add(
//         /*const*/ ISurface& surface,
//         /*const*/ IField& field,
//         /*const*/ SurfaceRegion& region,
//         int n_points,
//         const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix,
//         std::array<double, array_points> &xs, std::array<double, array_points> &ys, std::array<FieldValue, array_points> &FieldInRegion
//         )
// {
//     VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
//     VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());

//     // В тау-тау-эн;
//     VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
//     VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());



//     struct Task{
//             SurfaceRegion &region;
//             IField &field;
//             ISurface &surface;
//             VTKSurfaceSaver &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m;
//             int n_points;
//             std::array<double, array_points> &xs;
//             std::array<double, array_points> &ys;
//             std::array<FieldValue, array_points> &FieldInRegion;
//             size_t  i, j;
//             void operator()(){
//                 Vector2D xy(region.x_min + region.width() / (n_points-1) * i + region.width() / 2 / (n_points-1), region.y_min + region.height() / (n_points-1) * j + region.height() /2/(n_points-1));
//                 Position p = surface.point(xy);

//                 size_t count = i*(n_points-1) + j;

//                 xs[count]=xy[0];
//                 ys[count]=xy[1];

//                 // if (i==0) {std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;}
//                 // // std::cout << p[0] << " " << p[1] << " "<< p[2] << std::endl;

//                 FieldValue field_value = field.get(p);
//                 FieldInRegion[count]=field_value;
//                 //std::cout << "point ready" << std::endl;
//                 // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

//                 vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
//                 vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
//                 vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

//                 FieldValue field_rotated;
//                 // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
//                 // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//                 // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

//                 field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
//                 field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//                 field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

//                 // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
//                 // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));

//                 vtk_saver_surf_r.set_point(i, j, p, vec_modulus(field_value.E));
//                 vtk_saver_surf_m.set_point(i, j, p, vec_phases(field_value.E));
//             }
//             operator bool()const{
//                 return true;
//             }
//         };
//     using Threads = Threading::ThreadPool<ConcurrentQueue<Task>,14>;
//     typename Threads::TaskQueue tasks;
//     // using Threads = ThreadPool<std::function<void()>>;
//     // Threads::TaskQueue tasks;

//     for (int i = 0; i < n_points-1; i++)
//     {
//         for (int j = 0; j < n_points-1; j++)
//         {
//             tasks.push(Task{region, field, surface, vtk_saver_r, vtk_saver_i, vtk_saver_m, vtk_saver_surf_r, vtk_saver_surf_m, n_points, xs, ys, FieldInRegion, i, j});
//             // {
//             //     Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
//             //     Position p = surface.point(xy);
//             //     FieldValue field_value = field.get(p);
//             //     //std::cout << "point ready" << std::endl;
//             //     // field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

//             //     vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
//             //     vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
//             //     vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

//             //     FieldValue field_rotated;
//             //     // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
//             //     // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//             //     // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

//             //     field_rotated.E[0] = projection(field_value.E, surface.dS_over_dxdy(xy));
//             //     field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
//             //     field_rotated.E[2] = projection(field_value.E, surface.tau1(xy));

//             //     vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
//             //     vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
//             // });
//         }
//     }
//     Threads threads{tasks};
//     threads.run();
//     // Threads threads{tasks};
//     // threads.run();

//     vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
//     vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
//     vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

//     // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
//     // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());

//     vtk_saver_surf_r.save((filename_prefix + "_modulus_" + filename_suffix).c_str());
//     vtk_saver_surf_m.save((filename_prefix + "_phases_" + filename_suffix).c_str());
// }




// с добавкой фазы и старой параллелью
void plot_field_on_given_surface_mod(
        const ISurface& surface,
        const IField& field,
        const SurfaceRegion& region,
        double phase,
        int n_points,
        const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
{
    VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
    VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());

    // В тау-тау-эн;
    VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

    using Threads = ThreadPool<std::function<void()>>;
    Threads::TaskQueue tasks;

    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            tasks.push([&region, &field, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m, n_points, i, j, phase](){
                Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                Position p = surface.point(xy);
                FieldValue field_value = field.get(p);
                //std::cout << "point ready" << std::endl;
                field_value.E *=  exp(std::complex<double>(0.0, 1.0) * phase);

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
    Threads threads{tasks};
    threads.run();

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

    vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
}

//с учётом времени (ультракороткие) (старая параллель) (special-inverse-fourier.hpp)
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

    using Threads = ThreadPool<std::function<void()>>;
    Threads::TaskQueue tasks;

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
    Threads threads{tasks};
    threads.run();

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

    vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
}


//с учётом времени (ультракороткие) (3d array method - через fields-in-region) (старая параллель)
void plot_field_on_given_surface_with_time_3darr(
        const ISurface& surface,
        const SurfaceRegion& region,
        auto& fields, //FieldsInRegion<>& fields;
        std::vector<double> &freqs,
        std::vector<std::complex<double>> &amplF,
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

    fields.constructSummary(freqs, amplF, timing);

    using Threads = ThreadPool<std::function<void()>>;
    Threads::TaskQueue tasks;

    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            tasks.push([&region, &fields, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_surf_r, &vtk_saver_surf_m, n_points, timing, i, j](){
                Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                Position p = surface.point(xy);
                FieldValue field_value = fields.getSummary(i, j);
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
    Threads threads{tasks};
    threads.run();

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());

    vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
}

//3d array method (fields-in-region)
void plot_field_on_given_surface_3darr(
        const ISurface& surface,
        const SurfaceRegion& region,
        auto& fields, //FieldsInRegion<>& fields;
        // std::vector<double> &phases,
        // size_t iter,
        int n_points,
        const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
{
    VTKSurfaceSaver vtk_saver_r(n_points, n_points, (quantity_name + "_real").c_str());
    VTKSurfaceSaver vtk_saver_i(n_points, n_points, (quantity_name + "_imag").c_str());
    VTKSurfaceSaver vtk_saver_m(n_points, n_points, (quantity_name + "_max").c_str());
    VTKSurfaceSaver vtk_saver_int(n_points, n_points, (quantity_name + "_intens").c_str());

    // // В тау-тау-эн;
    // VTKSurfaceSaver vtk_saver_surf_r(n_points, n_points, (quantity_name + "_real").c_str());
    // VTKSurfaceSaver vtk_saver_surf_m(n_points, n_points, (quantity_name + "_max").c_str());

    using Threads = ThreadPool<std::function<void()>>;
    Threads::TaskQueue tasks;

    for (int i = 0; i < n_points; i++)
    {
        for (int j = 0; j < n_points; j++)
        {
            tasks.push([&region, &fields, &surface, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_int, /*&vtk_saver_surf_m,*/ n_points,  i, j/*, iter*/](){
                Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                Position p = surface.point(xy);
                // FieldValue field_value = fields.getHarmonic(0, i, j)
                FieldValue field_value = fields.getSummary(i, j);
                // FieldValue field_value = fields.getHarmonic(iter, i, j);
                //std::cout << "point ready "<< i << "," << j << std::endl;

                vtk_saver_r.set_point(i, j, p, vec_real(field_value.E));
                vtk_saver_i.set_point(i, j, p, vec_imag(field_value.E));
                vtk_saver_m.set_point(i, j, p, max_field(field_value.E));

                double absE;
                Vector ves;
                absE=0.0;
                ves=max_field(field_value.E);
                for (size_t l=0;l<3;l++){
                    absE+=ves[l]*ves[l];
                }
                absE=absE*4e13*50/16;
                ves[0]=absE;
                ves[1]=0;
                ves[2]=0;
                vtk_saver_int.set_point(i, j, p, ves);

                // FieldValue field_rotated;
                // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
                // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

                // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
                // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
            });
        }
    }
    Threads threads{tasks};
    threads.run();

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());
    vtk_saver_int.save((filename_prefix + "_intens_" + filename_suffix).c_str());

    // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
}

//3d-array method (fields-in-3dregion)
void plot_field_in_given_area_3darr(
        const IVolume& volume,
        const VolumeRegion& region,
        auto& fields, //FieldsInRegion<>& fields;
        // std::vector<double> &phases,
        // size_t iter,
        int n_points,
        const std::string& quantity_name, const std::string& filename_prefix, const std::string& filename_suffix)
{
    VTKVolumeSaver vtk_saver_r(n_points, n_points, n_points, (quantity_name + "_real").c_str());
    VTKVolumeSaver vtk_saver_i(n_points, n_points, n_points, (quantity_name + "_imag").c_str());
    VTKVolumeSaver vtk_saver_m(n_points, n_points, n_points, (quantity_name + "_max").c_str());
    VTKVolumeSaver vtk_saver_int(n_points, n_points, n_points, (quantity_name + "_intens").c_str());

    using Threads = ThreadPool<std::function<void()>>;
    Threads::TaskQueue tasks;

    for (int k = 0; k < n_points; k++)
    {
        for (int i = 0; i < n_points; i++)
        {
            for (int j = 0; j < n_points; j++)
            {
                tasks.push([&region, &fields, &volume, &vtk_saver_r, &vtk_saver_i, &vtk_saver_m, &vtk_saver_int, /*&vtk_saver_surf_m,*/ n_points,  i, j, k/*, iter*/](){
                    Vector xyz(region.x_min + region.widthx() / (n_points-1) * i, region.y_min + region.widthy() / (n_points-1) * j, region.z_min + region.height() / (n_points-1) * k);
                    Position p = volume.point(xyz);
                    FieldValue field_value = fields.getSummary(i, j, k);
                    // FieldValue field_value = fields.getHarmonic(iter, i, j);
                    //std::cout << "point ready "<< i << "," << j << std::endl;

                    vtk_saver_r.set_point(i, j, k, p, vec_real(field_value.E));
                    vtk_saver_i.set_point(i, j, k, p, vec_imag(field_value.E));
                    vtk_saver_m.set_point(i, j, k, p, max_field(field_value.E));

                    double absE;
                    Vector ves;
                    absE=0.0;
                    // ves=max_field(field_value.E);
                    ves=vec_modulus(field_value.E);
                    for (size_t l=0;l<3;l++){
                        absE+=ves[l]*ves[l];
                    }
                    absE=absE*4e13*50/16;
                    ves[0]=absE;
                    ves[1]=0;
                    ves[2]=0;
                    vtk_saver_int.set_point(i, j, k, p, ves);

                    // double absS, intens;
                    // Vector pointing_mod, intens_v, maxE, maxB, pointing;
                    // // VectorComplex pointing;
                    // absS=0.0;
                    // maxE=max_field(field_value.E);
                    // maxB=max_field(field_value.B);
                    // pointing=maxE%maxB;
                    // // pointing_mod=vec_modulus(pointing);
                    // for (size_t l=0;l<3;l++){
                    //     absS+=pointing[l]*pointing[l];
                    // }
                    // intens=sqrt(absS)*4e13*50/16;
                    // intens_v[0]=intens;
                    // intens_v[1]=0;
                    // intens_v[2]=0;
                    // vtk_saver_int.set_point(i, j, k, p, intens_v);

                    // FieldValue field_rotated;
                    // field_rotated.E[0] = projection(field_value.E, surface.tau1(xy));
                    // field_rotated.E[1] = projection(field_value.E, surface.tau2(xy));
                    // field_rotated.E[2] = projection(field_value.E, surface.dS_over_dxdy(xy));

                    // vtk_saver_surf_r.set_point(i, j, p, vec_real(field_rotated.E));
                    // vtk_saver_surf_m.set_point(i, j, p, max_field(field_rotated.E));
                });
            }
        }
    }
    Threads threads{tasks};
    threads.run();

    vtk_saver_r.save((filename_prefix + "_real_" + filename_suffix).c_str());
    vtk_saver_i.save((filename_prefix + "_imag_" + filename_suffix).c_str());
    vtk_saver_m.save((filename_prefix + "_max_" + filename_suffix).c_str());
    vtk_saver_int.save((filename_prefix + "_intens_" + filename_suffix).c_str());

    // vtk_saver_surf_r.save((filename_prefix + "_real_" + "surf_coords_" + filename_suffix).c_str());
    // vtk_saver_surf_m.save((filename_prefix + "_max_" + "surf_coords_" + filename_suffix).c_str());
}

#endif // VTKFIELDPLOTTER_HPP