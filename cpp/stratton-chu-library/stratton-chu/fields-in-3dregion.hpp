#ifndef FIELDSIN3DREGION_HPP
#define FIELDSIN3DREGION_HPP

#include "stratton-chu/types.hpp"
#include "stratton-chu/field.hpp"
#include "stratton-chu/surface.hpp"
#include "stratton-chu/volume.hpp"
#include "stratton-chu/utils.hpp"
#include <iostream>
#include "stratton-chu/stratton-chu-field.hpp"
#include "stratton-chu/spec-inverse-fourier.hpp"
#include <vector>
#include <list>
#include <deque>
#include <array>
#include <string_view>
#include <fstream>
#define IGNORE_BIND_CORE
#include "threading.hpp"

// constexpr double c = 29979245800.0;

template <const size_t n_points>
class FieldsIn3dRegion
{
public:
    FieldsIn3dRegion(IVolume& vol, VolumeRegion& region, std::string_view region_name, std::vector<StrattonChuReflectionLowPrec> &fields, const size_t n_harm):
        m_vol(vol), m_region(region), m_fields(fields), mn_harm(n_harm), mn_points(n_points), m_HarmonicsCatalogue(n_harm), m_SummField(n_points), m_SummFieldIntens(n_points), region_name(region_name)
    {   
        std::ifstream file_dump{this->region_name, std::ios::binary};
        if (file_dump.good())
        {
            // const std::size_t size = n_harm * sizeof(decltype(m_HarmonicsCatalogue)::value_type);
            for (auto& h :m_HarmonicsCatalogue)
            {
                file_dump.read(reinterpret_cast<char*>(h.data()), sizeof(typename decltype(m_HarmonicsCatalogue)::value_type));
            }
            // std::string field_file = "field_file.txt";
            // std::ofstream outr(field_file);
            // for (size_t w = 0; w < mn_harm; w++)
            // {
            //     for (size_t k = 0; k < mn_points; k++)
            //     {
            //         for (size_t i = 0; i < mn_points; i++)
            //         {
            //             for (size_t j = 0; j < mn_points; j++)
            //             {   
            //                 outr << "harm " << w << "   in point [" << i << ", " << j << ", " << k <<"] " << m_HarmonicsCatalogue[w][k][i][j].E[0] << " "<< m_HarmonicsCatalogue[w][k][i][j].E[1] << " "<< m_HarmonicsCatalogue[w][k][i][j].E[2] <<'\n';
            //             }
            //         }
            //     }
            // }          
            return;
        }
        file_dump.close();
        
        struct Task{
            FieldsIn3dRegion<n_points> *self;
            VolumeRegion &region;
            std::vector<StrattonChuReflectionLowPrec> &fields;
            IVolume &vol;
            size_t w, i, j, k;
            void operator()(){
                Vector xyz(region.x_min + region.widthx() / (n_points-1) * i, region.y_min + region.widthy() / (n_points-1) * j, region.z_min + region.height() / (n_points-1) * k);
                Position p = vol.point(xyz);
                // Position p = {0.0, 0.0, 0.0};
                self->m_HarmonicsCatalogue[w][k][i][j] = fields[w].get(p);
                // std::cout << self->m_HarmonicsCatalogue[w][k][i][j].E[0].real() << std::endl;
                // self->m_HarmonicsCatalogue[w][k][i][j].E[0].real(0.0);
                // std::cout << self->m_HarmonicsCatalogue[w][k][i][j].E[0].real() << std::endl;
                // std::cout << "point ready" << std::endl;
                if (i==0 && j==0 && k==0){
                    // std::cout << "point ready" << std::endl;
                    std::cout << "starting point " << w << std::endl;       
                }
                // if (i==80 && j==80){
                //     std::cout << "opa" << std::endl;      
                // } 
            }
            operator bool()const{
                return true;
            }
        };
        using Threads = Threading::ThreadPool<ConcurrentQueue<Task>,14>;
        typename Threads::TaskQueue tasks;
        // std::cout << "ready" << std::endl;

        if (n_harm==12){
            for (size_t w = 0; w < 2; w++)
            {
                for (size_t k = 0; k < n_points; k++)
                {
                    for (size_t i = 0; i < n_points; i++)
                    {
                        for (size_t j = 0; j < n_points; j++)
                        {   
                            tasks.push(Task{this, region, fields, vol, w, i, j, k});
                        }
                    }
                    // std::cout << "point "<<w<<" ready for " << region_name << std::endl;
                }
            }
        }
        else {
            for (size_t w = 0; w < n_harm; w++)
            {
                for (size_t k = 0; k < n_points; k++)
                {
                    for (size_t i = 0; i < n_points; i++)
                    {
                        for (size_t j = 0; j < n_points; j++)
                        {   
                            tasks.push(Task{this, region, fields, vol, w, i, j, k});
                        }
                    }
                    // std::cout << "point "<<w<<" ready for " << region_name << std::endl;
                }
            }
        }
        Threads threads{tasks};
        threads.run();
        if (n_harm==12){
            std::cout << "entering count simulation" << std::endl;

            for (size_t w = 2; w < n_harm; w++)
            {
                for (size_t k = 0; k < n_points; k++)
                {
                    for (size_t i = 0; i < n_points; i++)
                    {
                        for (size_t j = 0; j < n_points; j++)
                        {   
                            switch (w) {
                            case 2:   //1 otn y
                                // supv={1.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_HarmonicsCatalogue[w][k][i][j].E = m_HarmonicsCatalogue[1][k][i][n_points-1-j].E-2.0*(m_HarmonicsCatalogue[1][k][i][n_points-1-j].E*supv)*supv;                                
                                // m_HarmonicsCatalogue[w][k][i][j].B = m_HarmonicsCatalogue[1][k][i][n_points-1-j].B-2.0*(m_HarmonicsCatalogue[1][k][i][n_points-1-j].B*supv)*supv;   
                                
                                m_HarmonicsCatalogue[w][k][i][j].E[0].real(-m_HarmonicsCatalogue[1][k][n_points-1-i][j].E[0].real());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[0].imag(-m_HarmonicsCatalogue[1][k][n_points-1-i][j].E[0].imag());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[1].real(m_HarmonicsCatalogue[1][k][n_points-1-i][j].E[1].real());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[1].imag(m_HarmonicsCatalogue[1][k][n_points-1-i][j].E[1].imag());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[2].real(m_HarmonicsCatalogue[1][k][n_points-1-i][j].E[2].real());     //z
                                m_HarmonicsCatalogue[w][k][i][j].E[2].imag(m_HarmonicsCatalogue[1][k][n_points-1-i][j].E[2].imag());     //z
                                break;
                            case 3:   //0 otn y
                                // supv={1.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_HarmonicsCatalogue[w][k][i][j].E = m_HarmonicsCatalogue[0][k][i][n_points-1-j].E-2.0*(m_HarmonicsCatalogue[0][k][i][n_points-1-j].E*supv)*supv;                                
                                // m_HarmonicsCatalogue[w][k][i][j].B = m_HarmonicsCatalogue[0][k][i][n_points-1-j].B-2.0*(m_HarmonicsCatalogue[0][k][i][n_points-1-j].B*supv)*supv; 
                                m_HarmonicsCatalogue[w][k][i][j].E[0].real(-m_HarmonicsCatalogue[0][k][n_points-1-i][j].E[0].real());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[0].imag(-m_HarmonicsCatalogue[0][k][n_points-1-i][j].E[0].imag());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[1].real(m_HarmonicsCatalogue[0][k][n_points-1-i][j].E[1].real());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[1].imag(m_HarmonicsCatalogue[0][k][n_points-1-i][j].E[1].imag());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[2].real(m_HarmonicsCatalogue[0][k][n_points-1-i][j].E[2].real());     //z
                                m_HarmonicsCatalogue[w][k][i][j].E[2].imag(m_HarmonicsCatalogue[0][k][n_points-1-i][j].E[2].imag());     //z
                                break;
                            case 4:   //1 otn y,x
                                // supv={1.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_har_copy=m_HarmonicsCatalogue[1][k][n_points-1-i][n_points-1-j].E-2.0*(m_HarmonicsCatalogue[1][k][n_points-1-i][n_points-1-j].E*supv)*supv;
                                // supv={0.0+0.0*std::complex<double>(0.0, 1.0), 1.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_HarmonicsCatalogue[w][k][i][j].E = m_har_copy-2.0*(m_har_copy*supv)*supv;
                                // supv={1.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_har_copy=m_HarmonicsCatalogue[1][k][n_points-1-i][n_points-1-j].B-2.0*(m_HarmonicsCatalogue[1][k][n_points-1-i][n_points-1-j].B*supv)*supv;
                                // supv={0.0+0.0*std::complex<double>(0.0, 1.0), 1.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_HarmonicsCatalogue[w][k][i][j].B = m_har_copy-2.0*(m_har_copy*supv)*supv;
                                m_HarmonicsCatalogue[w][k][i][j].E[0].real(-m_HarmonicsCatalogue[1][k][n_points-1-i][n_points-1-j].E[0].real());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[0].imag(-m_HarmonicsCatalogue[1][k][n_points-1-i][n_points-1-j].E[0].imag());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[1].real(-m_HarmonicsCatalogue[1][k][n_points-1-i][n_points-1-j].E[1].real());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[1].imag(-m_HarmonicsCatalogue[1][k][n_points-1-i][n_points-1-j].E[1].imag());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[2].real(m_HarmonicsCatalogue[1][k][n_points-1-i][n_points-1-j].E[2].real());     //z
                                m_HarmonicsCatalogue[w][k][i][j].E[2].imag(m_HarmonicsCatalogue[1][k][n_points-1-i][n_points-1-j].E[2].imag());
                                break;
                            case 5:   //1 otn x
                                // supv={0.0+0.0*std::complex<double>(0.0, 1.0), 1.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_HarmonicsCatalogue[w][k][i][j].E = m_HarmonicsCatalogue[1][k][n_points-1-i][j].E-2.0*(m_HarmonicsCatalogue[1][k][n_points-1-i][j].E*supv)*supv;                                 
                                // m_HarmonicsCatalogue[w][k][i][j].B = m_HarmonicsCatalogue[1][k][n_points-1-i][j].B-2.0*(m_HarmonicsCatalogue[1][k][n_points-1-i][j].B*supv)*supv; 
                                m_HarmonicsCatalogue[w][k][i][j].E[0].real(m_HarmonicsCatalogue[1][k][i][n_points-1-j].E[0].real());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[0].imag(m_HarmonicsCatalogue[1][k][i][n_points-1-j].E[0].imag());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[1].real(-m_HarmonicsCatalogue[1][k][i][n_points-1-j].E[1].real());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[1].imag(-m_HarmonicsCatalogue[1][k][i][n_points-1-j].E[1].imag());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[2].real(m_HarmonicsCatalogue[1][k][i][n_points-1-j].E[2].real());     //z
                                m_HarmonicsCatalogue[w][k][i][j].E[2].imag(m_HarmonicsCatalogue[1][k][i][n_points-1-j].E[2].imag());
                                break;
                            case 6:   //0 otn y,z
                                // supv={0.0+0.0*std::complex<double>(0.0, 1.0), 1.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_HarmonicsCatalogue[w][k][i][j].E = 2.0*(m_HarmonicsCatalogue[0][n_points-1-k][i][n_points-1-j].E*supv)*supv-m_HarmonicsCatalogue[0][n_points-1-k][i][n_points-1-j].E;
                                // m_HarmonicsCatalogue[w][k][i][j].B = 2.0*(m_HarmonicsCatalogue[0][n_points-1-k][i][n_points-1-j].B*supv)*supv-m_HarmonicsCatalogue[0][n_points-1-k][i][n_points-1-j].B;
                                m_HarmonicsCatalogue[w][k][i][j].E[0].real(m_HarmonicsCatalogue[0][n_points-1-k][n_points-1-i][j].E[0].real());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[0].imag(m_HarmonicsCatalogue[0][n_points-1-k][n_points-1-i][j].E[0].imag());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[1].real(-m_HarmonicsCatalogue[0][n_points-1-k][n_points-1-i][j].E[1].real());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[1].imag(-m_HarmonicsCatalogue[0][n_points-1-k][n_points-1-i][j].E[1].imag());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[2].real(m_HarmonicsCatalogue[0][n_points-1-k][n_points-1-i][j].E[2].real());     //z
                                m_HarmonicsCatalogue[w][k][i][j].E[2].imag(m_HarmonicsCatalogue[0][n_points-1-k][n_points-1-i][j].E[2].imag());
                                break;
                            case 7:   //1 otn x,y,z
                                // supv={0.0+0.0*std::complex<double>(0.0, 1.0), 1.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_har_copy=2.0*(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][n_points-1-j].E*supv)*supv-m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][n_points-1-j].E;
                                // m_HarmonicsCatalogue[w][k][i][j].E = m_har_copy-2.0*(m_har_copy*supv)*supv;                                
                                // m_har_copy=2.0*(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][n_points-1-j].B*supv)*supv-m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][n_points-1-j].B;
                                // m_HarmonicsCatalogue[w][k][i][j].B = m_har_copy-2.0*(m_har_copy*supv)*supv;
                                m_HarmonicsCatalogue[w][k][i][j].E[0].real(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][n_points-1-j].E[0].real());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[0].imag(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][n_points-1-j].E[0].imag());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[1].real(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][n_points-1-j].E[1].real());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[1].imag(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][n_points-1-j].E[1].imag());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][n_points-1-j].E[2].real());     //z
                                m_HarmonicsCatalogue[w][k][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][n_points-1-j].E[2].imag());
                                break;
                            case 8:   //1 otn x,z
                                // supv={1.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_HarmonicsCatalogue[w][k][i][j].E = 2.0*(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][j].E*supv)*supv-m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][j].E;
                                // m_HarmonicsCatalogue[w][k][i][j].B = 2.0*(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][j].B*supv)*supv-m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][j].B;
                                m_HarmonicsCatalogue[w][k][i][j].E[0].real(-m_HarmonicsCatalogue[1][n_points-1-k][i][n_points-1-j].E[0].real());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[0].imag(-m_HarmonicsCatalogue[1][n_points-1-k][i][n_points-1-j].E[0].imag());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[1].real(m_HarmonicsCatalogue[1][n_points-1-k][i][n_points-1-j].E[1].real());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[1].imag(m_HarmonicsCatalogue[1][n_points-1-k][i][n_points-1-j].E[1].imag());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-k][i][n_points-1-j].E[2].real());     //z
                                m_HarmonicsCatalogue[w][k][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-k][i][n_points-1-j].E[2].imag());
                                break;
                            case 9:   //0 otn z
                                // supv={0.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0), 1.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_HarmonicsCatalogue[w][k][i][j].E = m_HarmonicsCatalogue[0][n_points-1-k][i][j].E-2.0*(m_HarmonicsCatalogue[0][n_points-1-k][i][j].E*supv)*supv;
                                // m_HarmonicsCatalogue[w][k][i][j].B = m_HarmonicsCatalogue[0][n_points-1-k][i][j].B-2.0*(m_HarmonicsCatalogue[0][n_points-1-k][i][j].B*supv)*supv;
                                m_HarmonicsCatalogue[w][k][i][j].E[0].real(-m_HarmonicsCatalogue[0][n_points-1-k][i][j].E[0].real());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[0].imag(-m_HarmonicsCatalogue[0][n_points-1-k][i][j].E[0].imag());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[1].real(-m_HarmonicsCatalogue[0][n_points-1-k][i][j].E[1].real());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[1].imag(-m_HarmonicsCatalogue[0][n_points-1-k][i][j].E[1].imag());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[2].real(m_HarmonicsCatalogue[0][n_points-1-k][i][j].E[2].real());     //z
                                m_HarmonicsCatalogue[w][k][i][j].E[2].imag(m_HarmonicsCatalogue[0][n_points-1-k][i][j].E[2].imag());
                                break;
                            case 10:   //1 otn z
                                // supv={0.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0), 1.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_HarmonicsCatalogue[w][k][i][j].E = m_HarmonicsCatalogue[1][n_points-1-k][i][j].E-2.0*(m_HarmonicsCatalogue[1][n_points-1-k][i][j].E*supv)*supv;
                                // m_HarmonicsCatalogue[w][k][i][j].B = m_HarmonicsCatalogue[1][n_points-1-k][i][j].B-2.0*(m_HarmonicsCatalogue[1][n_points-1-k][i][j].B*supv)*supv;
                                m_HarmonicsCatalogue[w][k][i][j].E[0].real(-m_HarmonicsCatalogue[1][n_points-1-k][i][j].E[0].real());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[0].imag(-m_HarmonicsCatalogue[1][n_points-1-k][i][j].E[0].imag());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[1].real(-m_HarmonicsCatalogue[1][n_points-1-k][i][j].E[1].real());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[1].imag(-m_HarmonicsCatalogue[1][n_points-1-k][i][j].E[1].imag());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-k][i][j].E[2].real());     //z
                                m_HarmonicsCatalogue[w][k][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-k][i][j].E[2].imag());
                                break;
                            case 11:   //1 otn y,z
                                // supv={0.0+0.0*std::complex<double>(0.0, 1.0), 1.0+0.0*std::complex<double>(0.0, 1.0), 0.0+0.0*std::complex<double>(0.0, 1.0)};
                                // m_HarmonicsCatalogue[w][k][i][j].E = 2.0*(m_HarmonicsCatalogue[1][n_points-1-k][i][n_points-1-j].E*supv)*supv-m_HarmonicsCatalogue[1][n_points-1-k][i][n_points-1-j].E;
                                // m_HarmonicsCatalogue[w][k][i][j].B = 2.0*(m_HarmonicsCatalogue[1][n_points-1-k][i][n_points-1-j].B*supv)*supv-m_HarmonicsCatalogue[1][n_points-1-k][i][n_points-1-j].B;
                                m_HarmonicsCatalogue[w][k][i][j].E[0].real(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][j].E[0].real());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[0].imag(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][j].E[0].imag());     //x
                                m_HarmonicsCatalogue[w][k][i][j].E[1].real(-m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][j].E[1].real());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[1].imag(-m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][j].E[1].imag());     //y
                                m_HarmonicsCatalogue[w][k][i][j].E[2].real(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][j].E[2].real());     //z
                                m_HarmonicsCatalogue[w][k][i][j].E[2].imag(m_HarmonicsCatalogue[1][n_points-1-k][n_points-1-i][j].E[2].imag());
                                break;
                            default:
                                break;                                
                            }
                            if (i==0 && j==0 && k==0){
                                std::cout << "starting point " << w << std::endl;       
                            }

                        }
                    }
                    // std::cout << "point "<<w<<" ready for " << region_name << std::endl;
                }
            }            
        }

        max_summary=0.0;
        irowx_max=0;
        irowy_max=0;
        icolumn_max=0;
        max_summary_exact=0.0;
        irowx_max_e=0;
        irowy_max_e=0;
        icolumn_max_e=0;
    }


    FieldValue getHarmonic(size_t ifreq, size_t irowx, size_t irowy, size_t icolumn) const 
    {
        return m_HarmonicsCatalogue[ifreq][icolumn][irowx][irowy];
    }


    void constructSummary(std::vector<double> &freqs, std::vector<std::complex<double>> &amplF, const double time)
    {   

        // std::cout << "freqs "<< freqs.size() << std::endl;

        // std::cout << "ampls "<< amplF.size() << std::endl;
        // std::cout << "mharm "<< m_HarmonicsCatalogue.size() << std::endl;
        mono_flag=false;
        FieldValue null_field;
        for (size_t k = 0; k < mn_points; k++)
        {
            for (size_t i = 0; i < mn_points; i++)
            {
                for (size_t j = 0; j < mn_points; j++)
                {
                    m_SummField[k][i][j]=null_field;
                    for (size_t w = 0; w < mn_harm; w++)
                    {                       
                        //std::cout << amplF[w] << " "<< freqs[w] << std::endl;
                        m_SummField[k][i][j].E += m_HarmonicsCatalogue[w][k][i][j].E*amplF[w]*exp(std::complex<double>(0.0, 1.0) * freqs[w] * time);                        
                    }
                    // std::cout << "hihihi" << std::endl;
                }           
            }
        }

    }

    void constructSummaryMono(std::vector<double> &phases)
    {   
        // std::cout << "freqs "<< freqs.size() << std::endl;
        // std::cout << "ampls "<< amplF.size() << std::endl;
        // std::cout << "mharm "<< m_HarmonicsCatalogue.size() << std::endl;

        mono_flag=true;
        FieldValue null_field;
        for (size_t k = 0; k < mn_points; k++)
        {
            for (size_t i = 0; i < mn_points; i++)
            {
                for (size_t j = 0; j < mn_points; j++)
                {
                    // m_SummField[k][i][j]=null_field;
                    m_SummField[k][i][j].E = m_HarmonicsCatalogue[0][k][i][j].E*exp(std::complex<double>(0.0, 1.0) * phases[0]);
                    m_SummField[k][i][j].B = m_HarmonicsCatalogue[0][k][i][j].B*exp(std::complex<double>(0.0, 1.0) * phases[0]);
                    // if (mn_harm>0) {
                    if (mn_harm>1) {
                        for (size_t w = 1; w < mn_harm; w++)
                        {                       
                            //std::cout << amplF[w] << " "<< freqs[w] << std::endl;
                            m_SummField[k][i][j].E += m_HarmonicsCatalogue[w][k][i][j].E*exp(std::complex<double>(0.0, 1.0) * phases[w]);
                            m_SummField[k][i][j].B += m_HarmonicsCatalogue[w][k][i][j].B*exp(std::complex<double>(0.0, 1.0) * phases[w]);                        
                        }
                        // std::cout << "hihihi" << std::endl;
                    }
                }           
            }
        }
    }

    FieldValue getSummary(size_t irowx, size_t irowy, size_t icolumn) const 
    {
        return m_SummField[icolumn][irowx][irowy];
    }

    double getLambda() const 
    {
        if (mono_flag) {return m_fields[0].lambda();}
        else {
            std::cout << "non-monochromatic field, particular lambda is unknown" << std::endl;
            return 0.0;
        }
        
    }

    void findSummaryMax(){
        double absE,maxabsE;
        size_t irx_max, iry_max, ic_max;
        Vector ves;
        irx_max=0;
        iry_max=0;
        ic_max=0;
        maxabsE=0.0;
        for (size_t k = 0; k < mn_points; k++)
        {
            for (size_t i = 0; i < mn_points; i++)
            {
                for (size_t j = 0; j < mn_points; j++)
                {
                    absE=0.0;
                    ves=vec_modulus(m_SummField[k][i][j].E);
                    for (size_t l=0;l<3;l++){
                        absE+=ves[l]*ves[l];
                    }
                    absE=sqrt(absE);
                    // std::cout << absE << std::endl;
                    if (absE>maxabsE){
                        maxabsE=absE;
                        irx_max=i;
                        iry_max=j;
                        ic_max=k;
                        // std::cout << "q" << std::endl;
                    }
                }           
            }
        }
        max_summary=maxabsE;
        irowx_max=irx_max;
        irowy_max=iry_max;
        icolumn_max=ic_max;
    }

    void findSummaryMaxExact(){
        double absE,maxabsE;
        size_t irx_max, iry_max, ic_max;
        Vector ves;
        irx_max=0;
        iry_max=0;
        ic_max=0;
        maxabsE=0.0;
        for (size_t k = 0; k < mn_points; k++)
        {
            for (size_t i = 0; i < mn_points; i++)
            {
                for (size_t j = 0; j < mn_points; j++)
                {
                    absE=0.0;
                    ves=max_field(m_SummField[k][i][j].E);
                    for (size_t l=0;l<3;l++){
                        absE+=ves[l]*ves[l];
                    }
                    absE=sqrt(absE);
                    // std::cout << absE << std::endl;
                    if (absE>maxabsE){
                        maxabsE=absE;
                        irx_max=i;
                        iry_max=j;
                        ic_max=k;
                        // std::cout << "q" << std::endl;
                    }
                }           
            }
        }
        // ves=vec_modulus(m_SummField[150][150].E);
        // for (size_t l;l<3;l++){
        //     std::cout << ves[l] << std::endl;
        // }
        max_summary_exact=maxabsE;
        irowx_max_e=irx_max;
        irowy_max_e=iry_max;
        icolumn_max_e=ic_max;

    }

    bool checkMax(){
        if (irowx_max_e==irowx_max && irowy_max_e==irowy_max && icolumn_max_e==icolumn_max){
            return true;
        } else {
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

    Vector MaxIntensityPoint()
    {   
        Vector xyz(m_region.x_min + m_region.widthx() / (mn_points-1) * irowx_max, m_region.y_min + m_region.widthy() / (mn_points-1) * irowy_max, m_region.z_min + m_region.height() / (mn_points-1) * icolumn_max);
        return xyz;
    }


    ~FieldsIn3dRegion()
    {
        std::fstream file_dump(region_name, std::ios::out | std::ios::binary);
        for (auto& h :m_HarmonicsCatalogue)
        {
            file_dump.write(reinterpret_cast<const char*>(h.data()), sizeof(typename decltype(m_HarmonicsCatalogue)::value_type));
        }            
    }

    

private:
    IVolume& m_vol;
    VolumeRegion m_region;
    std::vector<StrattonChuReflectionLowPrec> &m_fields;    
    size_t mn_harm;
    size_t mn_points;
    double max_summary;
    size_t irowx_max, irowy_max, icolumn_max;
    double max_summary_exact;
    size_t irowx_max_e, irowy_max_e, icolumn_max_e;
    //std::array<std::vector<std::vector<FieldValue>>, harmonics_number> m_HarmonicsCatalogue (mn_points, std::vector<FieldValue>(mn_points));

    //std::vector<std::vector<FieldValue> m_HarmonicsCatalogue (n_points, std::vector<FieldValue>(n_points))[harmonics_number];
    
    // std::vector<std::array<std::array<std::array<FieldValue, n_points>, n_points>, n_points>> m_HarmonicsCatalogue;
    // std::vector<std::array<std::array<FieldValue, n_points>, n_points>> m_HarmonicsCatalogue;
    std::deque<std::array<std::array<std::array<FieldValue, n_points>, n_points>, n_points>> m_HarmonicsCatalogue;

    std::vector<std::array<std::array<FieldValue, n_points>, n_points>> m_SummField;
    std::vector<std::array<double, n_points>> m_SummFieldIntens;

    std::string region_name;
    
    bool mono_flag;
};



#endif // FIELDSIN3DREGION_HPP