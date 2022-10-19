#include "stratton-chu/field-in-region.hpp"

#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <array>
#include "constants.hpp"
#include "thread_pool.hpp"

FieldInRegion::FieldInRegion(ISurface& surf, SurfaceRegion& region, std::vector<StrattonChuReflection> &fields, int n_harmonics, int n_points):
    m_surf(surf), m_region(region), m_fields(fields), mn_harm(n_harmonics), mn_points(n_points)
{
    std::vector<std::vector<std::vector<FieldValue>>> FieldCatalogue(n_harmonics, std::vector<std::vector<FieldValue>>(n_points, std::vector<FieldValue>(n_points)));
    function task{[&FieldCatalogue, &region, &fields, &surf, n_points, w = 0, i = 0, j = 0](){}};
    using TaskQueue = ConcurentQueue<decltype(task)>;
    using ThreadPool = ThreadPoolImpl<TaskQueue>;
    
    TaskQueue tasks;
    for (int w = 0; w<n_harmonics; w++)
    {
        for (int i = 0; i < n_points; i++)
        {
            for (int j = 0; j < n_points; j++)
            {   
                tasks.push([&FieldCatalogue, &region, &fields, &surf, n_points, w, i, j](){
                    Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                    Position p = surf.point(xy);
                    FieldCatalogue[w][i][j] = fields[w].get(p);
                    //std::cout << "point ready" << std::endl;
                });
            }
        }
    }
    ThreadPool threads{tasks};
    threads.run(worker<TaskQueue::value_type>);
    
}

FieldValue FieldInRegion::get(int ifreq, int irow, int icolumn) const 
{
    return m_FieldCatalogue[ifreq][irow][icolumn];
}

std::vector<std::vector<FieldValue>> FieldInRegion::getIFT(std::vector<double> &freqs, std::vector<std::complex<double>> &amplF, const double time)
{
    std::vector<std::vector<FieldValue>> result(mn_points, std::vector<FieldValue>(mn_points)); 

    for (int i = 0; i < mn_points; i++)
    {
        for (int j = 0; j < mn_points; j++)
        {
            for (int w = 0; w < mn_harm; w++)
            {
                    Vector2D xy(m_region.x_min + m_region.width() / (mn_points-1) * i, m_region.y_min + m_region.height() / (mn_points-1) * j);
                    Position p = m_surf.point(xy);
                    result[i][j].E=result[i][j].E+m_FieldCatalogue[w][i][j].E*amplF[w]*exp(std::complex<double>(0.0, 1.0) * freqs[w] * time);
                    //std::cout << "point ready" << std::endl;                
            }

        }
    }

    return result;
}