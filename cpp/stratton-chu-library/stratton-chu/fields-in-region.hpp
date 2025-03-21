#ifndef FIELDSINREGION_HPP
#define FIELDSINREGION_HPP

#include "stratton-chu/types.hpp"
#include "stratton-chu/field.hpp"
#include "stratton-chu/surface.hpp"
#include "stratton-chu/utils.hpp"
#include <iostream>
#include "stratton-chu/stratton-chu-field.hpp"
#include "stratton-chu/spec-inverse-fourier.hpp"
#include <vector>
#include <array>
#include <list>
#include <deque>
#include <string_view>
#include <fstream>
//#include "thread_pool.hpp"
#define IGNORE_BIND_CORE
#include "threading.hpp"

// constexpr double c = 29979245800.0;

template <const size_t n_points>
class FieldsInRegion
{
public:
    FieldsInRegion(ISurface& surf, SurfaceRegion& region, std::string_view region_name, std::vector<StrattonChuReflectionLowPrec> &fields, const size_t n_harm):
        m_surf(surf), m_region(region), m_fields(fields), mn_harm(n_harm), mn_points(n_points), m_HarmonicsCatalogue(n_harm), m_SummField(n_points), m_SummFieldIntens(n_points), region_name(region_name)
    {   
        std::ifstream file_dump{this->region_name, std::ios::binary};
        if (file_dump.good())
        {
            for (auto& h :m_HarmonicsCatalogue)
            {
                file_dump.read(reinterpret_cast<char*>(h.data()), sizeof(typename decltype(m_HarmonicsCatalogue)::value_type));
            }            
            return;
            // const std::size_t size = n_harm * sizeof(m_HarmonicsCatalogue[0]);
            // file_dump.read(reinterpret_cast<char*>(m_HarmonicsCatalogue.data()), size);
            // return;
        }
        file_dump.close();

        struct Task{
            FieldsInRegion<n_points> *self;
            SurfaceRegion &region;
            std::vector<StrattonChuReflectionLowPrec> &fields;
            ISurface &surf;
            size_t w, i, j;
            void operator()(){
                Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                Position p = surf.point(xy);
                self->m_HarmonicsCatalogue[w][i][j] = fields[w].get(p);
                // if (i==0 && j==0){
                //     std::cout << "start point ready " << w << std::endl;      
                // }
                // if (i==n_points-1 && j==n_points-1){
                //     std::cout << "end point ready" << std::endl;      
                // }
                // if (i==j){
                //     std::cout << "point ready " << i << std::endl;      
                // }
                // std::cout << "point ready " << std::endl;
                std::cout << "point " << i*n_points+j << " ready" << std::endl;                 
            }
            operator bool()const{
                return true;
            }
        };
        using Threads = Threading::ThreadPool<ConcurrentQueue<Task>,14>;
        // using Threads = ThreadPool<Task>;
        typename Threads::TaskQueue tasks;

        // for (size_t j = 0; j < n_points; j++)
        // {   
        //     tasks.push(Task{this, region, fields, surf, 0, 0, j});
        //     tasks.push(Task{this, region, fields, surf, 0, n_points-1, j});
        //     if (j!=0 && j!=n_points-1) {
        //         tasks.push(Task{this, region, fields, surf, 0, j, 0});
        //         tasks.push(Task{this, region, fields, surf, 0, j, n_points-1});
        //     }
        // }

        // for (size_t j = 1; j < n_points-1; j++)
        // {   
        //     tasks.push(Task{this, region, fields, surf, 0, 1, j});
        //     tasks.push(Task{this, region, fields, surf, 0, n_points-2, j});
        //     if (j!=1 && j!=n_points-2) {
        //         tasks.push(Task{this, region, fields, surf, 0, j, 1});
        //         tasks.push(Task{this, region, fields, surf, 0, j, n_points-2});
        //     }
        // }

        // for (size_t j = 2; j < n_points-2; j++)
        // {   
        //     tasks.push(Task{this, region, fields, surf, 0, 2, j});
        //     tasks.push(Task{this, region, fields, surf, 0, n_points-3, j});
        //     if (j!=2 && j!=n_points-3) {
        //         tasks.push(Task{this, region, fields, surf, 0, j, 2});
        //         tasks.push(Task{this, region, fields, surf, 0, j, n_points-3});
        //     }
        // }

        // for (size_t i = 3; i < n_points-3; i++)
        // {
        //     for (size_t j = 3; j < n_points-3; j++)
        //     {   
        //         tasks.push(Task{this, region, fields, surf, 0, i, j});
        //     }
        // }


        // for (size_t w = 0; w < n_harm; w++)
        // {
        //     for (size_t i = 0; i < n_points; i++)
        //     {
        //         for (size_t j = 0; j < n_points; j++)
        //         {   
        //             tasks.push(Task{this, region, fields, surf, w, i, j});
        //         }
        //     }
        //     // std::cout << "point "<<w<<" ready for " << region_name << std::endl;
        // }


        if (n_harm==12){
            // for (size_t w = 0; w < 2; w++)
            for (size_t w = 0; w < mn_harm; w++)
            {   
                for (size_t j = 0; j < n_points; j++)
                {   
                    tasks.push(Task{this, region, fields, surf, w, 0, j});
                    tasks.push(Task{this, region, fields, surf, w, n_points-1, j});
                    if (j!=0 && j!=n_points-1) {
                        tasks.push(Task{this, region, fields, surf, w, j, 0});
                        tasks.push(Task{this, region, fields, surf, w, j, n_points-1});
                    }
                }

                for (size_t j = 1; j < n_points-1; j++)
                {   
                    tasks.push(Task{this, region, fields, surf, w, 1, j});
                    tasks.push(Task{this, region, fields, surf, w, n_points-2, j});
                    if (j!=1 && j!=n_points-2) {
                        tasks.push(Task{this, region, fields, surf, w, j, 1});
                        tasks.push(Task{this, region, fields, surf, w, j, n_points-2});
                    }
                }

                for (size_t j = 2; j < n_points-2; j++)
                {   
                    tasks.push(Task{this, region, fields, surf, w, 2, j});
                    tasks.push(Task{this, region, fields, surf, w, n_points-3, j});
                    if (j!=2 && j!=n_points-3) {
                        tasks.push(Task{this, region, fields, surf, w, j, 2});
                        tasks.push(Task{this, region, fields, surf, w, j, n_points-3});
                    }
                }

                for (size_t i = 3; i < n_points-3; i++)
                {
                    for (size_t j = 3; j < n_points-3; j++)
                    {   
                        tasks.push(Task{this, region, fields, surf, w, i, j});
                    }
                }
                // for (size_t i = 0; i < n_points; i++)
                // {
                //     for (size_t j = 0; j < n_points; j++)
                //     {   
                //         tasks.push(Task{this, region, fields, surf, w, i, j});
                //     }
                // }
            }
        }
        else {
            for (size_t w = 0; w < n_harm; w++)
            {
                for (size_t i = 0; i < n_points; i++)
                {
                    for (size_t j = 0; j < n_points; j++)
                    {   
                        tasks.push(Task{this, region, fields, surf, w, i, j});
                    }
                }
            }
        }
        Threads threads{tasks};
        threads.run();

        // if (n_harm==12){
        //     std::cout << "entering count simulation" << std::endl;

        //     for (size_t w = 2; w < n_harm; w++)
        //     {
        //         if (region_name=="focal_plane_long_12-beams_new_config")
        //         {
        //             for (size_t i = 0; i < n_points; i++)
        //             {
        //                 for (size_t j = 0; j < n_points; j++)
        //                 {   
        //                     switch (w) {
        //                     case 2:   //1 otn y
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[1][n_points-1-i][j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[1][n_points-1-i][j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(m_HarmonicsCatalogue[1][n_points-1-i][j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(m_HarmonicsCatalogue[1][n_points-1-i][j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-i][j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-i][j].E[2].imag());     //z
        //                         break;
        //                     case 3:   //0 otn y
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[0][n_points-1-i][j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[0][n_points-1-i][j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(m_HarmonicsCatalogue[0][n_points-1-i][j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(m_HarmonicsCatalogue[0][n_points-1-i][j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[0][n_points-1-i][j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[0][n_points-1-i][j].E[2].imag());     //z
        //                         break;
        //                     case 4:   //1 otn y,x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[2].imag());
        //                         break;
        //                     case 5:   //1 otn x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(m_HarmonicsCatalogue[1][i][n_points-1-j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(m_HarmonicsCatalogue[1][i][n_points-1-j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[1][i][n_points-1-j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[1][i][n_points-1-j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][i][n_points-1-j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][i][n_points-1-j].E[2].imag());
        //                         break;
        //                     case 6:   //0 otn y,z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(m_HarmonicsCatalogue[0][n_points-1-i][j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(m_HarmonicsCatalogue[0][n_points-1-i][j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[0][n_points-1-i][j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[0][n_points-1-i][j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[0][n_points-1-i][j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[0][n_points-1-i][j].E[2].imag());
        //                         break;
        //                     case 7:   //1 otn x,y,z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[2].imag());
        //                         break;
        //                     case 8:   //1 otn x,z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[1][i][n_points-1-j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[1][i][n_points-1-j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(m_HarmonicsCatalogue[1][i][n_points-1-j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(m_HarmonicsCatalogue[1][i][n_points-1-j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][i][n_points-1-j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][i][n_points-1-j].E[2].imag());
        //                         break;
        //                     case 9:   //0 otn z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[0][i][j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[0][i][j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[0][i][j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[0][i][j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[0][i][j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[0][i][j].E[2].imag());
        //                         break;
        //                     case 10:   //1 otn z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[1][i][j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[1][i][j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[1][i][j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[1][i][j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][i][j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][i][j].E[2].imag());
        //                         break;
        //                     case 11:   //1 otn y,z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(m_HarmonicsCatalogue[1][n_points-1-i][j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(m_HarmonicsCatalogue[1][n_points-1-i][j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[1][n_points-1-i][j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[1][n_points-1-i][j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-i][j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-i][j].E[2].imag());
        //                         break;
        //                     default:
        //                         break;                                
        //                     }
        //                     if (i==0 && j==0){
        //                         std::cout << "starting point " << w << std::endl;       
        //                     }

        //                 }
        //             }
        //             // std::cout << "point "<<w<<" ready for " << region_name << std::endl;
        //         }
        //         if (region_name=="focal_plane_trans_12-beams_new_config")
        //             for (size_t i = 0; i < n_points; i++)
        //             {
        //                 for (size_t j = 0; j < n_points; j++)
        //                 {   
        //                     switch (w) {
        //                     case 2:   //1 otn y
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[1][n_points-1-i][j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[1][n_points-1-i][j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(m_HarmonicsCatalogue[1][n_points-1-i][j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(m_HarmonicsCatalogue[1][n_points-1-i][j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-i][j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-i][j].E[2].imag());     //z
        //                         break;
        //                     case 3:   //0 otn y
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[0][n_points-1-i][j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[0][n_points-1-i][j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(m_HarmonicsCatalogue[0][n_points-1-i][j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(m_HarmonicsCatalogue[0][n_points-1-i][j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[0][n_points-1-i][j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[0][n_points-1-i][j].E[2].imag());     //z
        //                         break;
        //                     case 4:   //1 otn y,x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[1][n_points-1-i][j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[1][n_points-1-i][j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[1][n_points-1-i][j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[1][n_points-1-i][j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-i][j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-i][j].E[2].imag());
        //                         break;
        //                     case 5:   //1 otn x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(m_HarmonicsCatalogue[1][i][j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(m_HarmonicsCatalogue[1][i][j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[1][i][j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[1][i][j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][i][j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][i][j].E[2].imag());
        //                         break;
        //                     case 6:   //0 otn y,z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(m_HarmonicsCatalogue[0][n_points-1-i][n_points-1-j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(m_HarmonicsCatalogue[0][n_points-1-i][n_points-1-j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[0][n_points-1-i][n_points-1-j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[0][n_points-1-i][n_points-1-j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[0][n_points-1-i][n_points-1-j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[0][n_points-1-i][n_points-1-j].E[2].imag());
        //                         break;
        //                     case 7:   //1 otn x,y,z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[2].imag());
        //                         break;
        //                     case 8:   //1 otn x,z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[1][i][n_points-1-j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[1][i][n_points-1-j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(m_HarmonicsCatalogue[1][i][n_points-1-j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(m_HarmonicsCatalogue[1][i][n_points-1-j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][i][n_points-1-j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][i][n_points-1-j].E[2].imag());
        //                         break;
        //                     case 9:   //0 otn z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[0][i][n_points-1-j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[0][i][n_points-1-j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[0][i][n_points-1-j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[0][i][n_points-1-j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[0][i][n_points-1-j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[0][i][n_points-1-j].E[2].imag());
        //                         break;
        //                     case 10:   //1 otn z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(-m_HarmonicsCatalogue[1][i][n_points-1-j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(-m_HarmonicsCatalogue[1][i][n_points-1-j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[1][i][n_points-1-j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[1][i][n_points-1-j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][i][n_points-1-j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][i][n_points-1-j].E[2].imag());
        //                         break;
        //                     case 11:   //1 otn y,z
        //                         m_HarmonicsCatalogue[w][i][j].E[0].real(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[0].real());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[0].imag(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[0].imag());     //x
        //                         m_HarmonicsCatalogue[w][i][j].E[1].real(-m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[1].real());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[1].imag(-m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[1].imag());     //y
        //                         m_HarmonicsCatalogue[w][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[2].real());     //z
        //                         m_HarmonicsCatalogue[w][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-i][n_points-1-j].E[2].imag());
        //                         break;
        //                     default:
        //                         break;                                
        //                     }
        //                     if (i==0 && j==0){
        //                         std::cout << "starting point " << w << std::endl;       
        //                     }

        //                 }
        //             }
        //         }
        //     }            
        // }
        max_summary=0.0;
        irow_max=0;
        icolumn_max=0;
        max_summary_exact=0.0;
        irow_max_e=0;
        icolumn_max_e=0;
    }


    FieldValue getHarmonic(size_t ifreq, size_t irow, size_t icolumn) const 
    {
        return m_HarmonicsCatalogue[ifreq][irow][icolumn];
    }


    void constructSummary(std::vector<double> &freqs, std::vector<std::complex<double>> &amplF, const double time)
    {   

        // std::cout << "freqs "<< freqs.size() << std::endl;

        // std::cout << "ampls "<< amplF.size() << std::endl;
        // std::cout << "mharm "<< m_HarmonicsCatalogue.size() << std::endl;
        
        FieldValue null_field;
        for (size_t i = 0; i < mn_points; i++)
        {
            for (size_t j = 0; j < mn_points; j++)
            {
                m_SummField[i][j]=null_field;
                for (size_t w = 0; w < mn_harm; w++)
                {                       
                    //std::cout << amplF[w] << " "<< freqs[w] << std::endl;
                    m_SummField[i][j].E += m_HarmonicsCatalogue[w][i][j].E*amplF[w]*exp(std::complex<double>(0.0, 1.0) * freqs[w] * time);                        
                }
                // std::cout << "hihihi" << std::endl;
            }           
        }

    }

    void constructSummaryMono(std::vector<double> &phases)
    {   
        // std::cout << "freqs "<< freqs.size() << std::endl;
        // std::cout << "ampls "<< amplF.size() << std::endl;
        // std::cout << "mharm "<< m_HarmonicsCatalogue.size() << std::endl;
        
        FieldValue null_field;
        for (size_t i = 0; i < mn_points; i++)
        {
            for (size_t j = 0; j < mn_points; j++)
            {
                m_SummField[i][j]=null_field;
                for (size_t w = 0; w < mn_harm; w++)
                {                       
                    //std::cout << amplF[w] << " "<< freqs[w] << std::endl;
                    m_SummField[i][j].E += m_HarmonicsCatalogue[w][i][j].E*exp(std::complex<double>(0.0, 1.0) * phases[w]);
                    m_SummField[i][j].B += m_HarmonicsCatalogue[w][i][j].B*exp(std::complex<double>(0.0, 1.0) * phases[w]);                        
                }
                // std::cout << "hihihi" << std::endl;
            }           
        }
    }

    FieldValue getSummary(size_t irow, size_t icolumn) const 
    {
        return m_SummField[irow][icolumn];
    }

    void findSummaryMax(){
        double absE,maxabsE;
        size_t ir_max, ic_max;
        Vector ves;
        ir_max=0;
        ic_max=0;
        maxabsE=0.0;
        for (size_t i = 0; i < mn_points; i++)
        {
            for (size_t j = 0; j < mn_points; j++)
            {
                absE=0.0;
                ves=vec_modulus(m_SummField[i][j].E);
                for (size_t l=0;l<3;l++){
                    absE+=ves[l]*ves[l];
                }
                absE=sqrt(absE);
                // std::cout << absE << std::endl;
                if (absE>maxabsE){
                    maxabsE=absE;
                    ir_max=i;
                    ic_max=j;
                    // std::cout << "q" << std::endl;
                }
            }           
        }
        max_summary=maxabsE;
        irow_max=ir_max;
        icolumn_max=ic_max;
    }

    void findSummaryMaxExact(){
        double absE,maxabsE;
        size_t ir_max, ic_max;
        Vector ves;
        ir_max=0;
        ic_max=0;
        maxabsE=0.0;
        for (size_t i = 0; i < mn_points; i++)
        {
            for (size_t j = 0; j < mn_points; j++)
            {
                absE=0.0;
                ves=max_field(m_SummField[i][j].E);
                for (size_t l=0;l<3;l++){
                    absE+=ves[l]*ves[l];
                }
                absE=sqrt(absE);
                // std::cout << absE << std::endl;
                if (absE>maxabsE){
                    maxabsE=absE;
                    ir_max=i;
                    ic_max=j;
                    // std::cout << "q" << std::endl;
                }
            }           
        }
        // ves=vec_modulus(m_SummField[150][150].E);
        // for (size_t l;l<3;l++){
        //     std::cout << ves[l] << std::endl;
        // }
        max_summary_exact=maxabsE;
        irow_max_e=ir_max;
        icolumn_max_e=ic_max;

    }

    bool checkMax(){
        if (irow_max_e==irow_max && icolumn_max_e==icolumn_max){
            return true;
        }
        else{
            return false;
        }
    }

    double calculateMaxIntensitySummary()
    {   
        // double c = 29979245800.0;
        double intensity;
        // intensity = c/8/M_PI*max_summary*max_summary;
        intensity = max_summary*max_summary*4e13*50/16;
        return intensity;
    }

    double calculateTrueMaxIntensitySummary()
    {   
        // double c = 29979245800.0;
        double intensity;
        // intensity = c/8/M_PI*max_summary*max_summary;
        intensity = max_summary_exact*max_summary_exact*4e13*50/16;
        return intensity;
    }

    Vector2D MaxIntensityPoint()
    {   
        Vector2D xy(m_region.x_min + m_region.width() / (mn_points-1) * irow_max, m_region.y_min + m_region.width() / (mn_points-1) * icolumn_max);
        return xy;
    }



    ~FieldsInRegion()
    {
        const std::size_t size = m_HarmonicsCatalogue.size() * sizeof(m_HarmonicsCatalogue[0]);
        std::fstream file_dump(region_name, std::ios::out | std::ios::binary);
        file_dump.write(reinterpret_cast<const char*>(m_HarmonicsCatalogue.data()), size);
    }

    

private:
    ISurface& m_surf;
    SurfaceRegion m_region;
    std::vector<StrattonChuReflectionLowPrec> &m_fields;    
    size_t mn_harm;
    size_t mn_points;
    double max_summary;
    size_t irow_max, icolumn_max;
    double max_summary_exact;
    size_t irow_max_e, icolumn_max_e;
    //std::array<std::vector<std::vector<FieldValue>>, harmonics_number> m_HarmonicsCatalogue (mn_points, std::vector<FieldValue>(mn_points));

    //std::vector<std::vector<FieldValue> m_HarmonicsCatalogue (n_points, std::vector<FieldValue>(n_points))[harmonics_number];
    
    std::vector<std::array<std::array<FieldValue, n_points>, n_points>> m_HarmonicsCatalogue;

    std::vector<std::array<FieldValue, n_points>> m_SummField;
    std::vector<std::array<double, n_points>> m_SummFieldIntens;

    std::string region_name;
};



#endif // FIELDSINREGION_HPP