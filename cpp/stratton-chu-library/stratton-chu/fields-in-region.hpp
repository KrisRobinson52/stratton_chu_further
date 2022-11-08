#ifndef FIELDSINREGION_HPP
#define FIELDSINREGION_HPP

#include "stratton-chu/types.hpp"
#include "stratton-chu/field.hpp"
#include "stratton-chu/surface.hpp"
#include <iostream>
#include "stratton-chu/stratton-chu-field.hpp"
#include "stratton-chu/spec-inverse-fourier.hpp"
#include <vector>
#include <array>
#include "thread_pool.hpp"

template <const size_t n_points>
class FieldsInRegion
{
public:
    FieldsInRegion(ISurface& surf, SurfaceRegion& region, std::vector<StrattonChuReflection> &fields, const size_t n_harm):
        m_surf(surf), m_region(region), m_fields(fields), mn_harm(n_harm), m_HarmonicsCatalogue(n_harm),m_SummField(n_points)
    {   
        //std::vector<std::vector<FieldValue> m_HarmonicsCatalogue (n_points, std::vector<FieldValue>(n_points))[harmonics_number];


        function task{[this, &region, &fields, &surf, w = 0lu, i = 0lu, j = 0lu](){}};
        using TaskQueue = ConcurentQueue<decltype(task)>;
        using ThreadPool = ThreadPoolImpl<TaskQueue>;
        
        TaskQueue tasks;
        for (size_t w = 0; w < n_harm; w++)
        {
            for (size_t i = 0; i < n_points; i++)
            {
                for (size_t j = 0; j < n_points; j++)
                {   
                    tasks.push([this, &region, &fields, &surf, w, i, j](){
                        Vector2D xy(region.x_min + region.width() / (n_points-1) * i, region.y_min + region.height() / (n_points-1) * j);
                        Position p = surf.point(xy);
                        m_HarmonicsCatalogue[w][i][j] = fields[w].get(p);
                        //std::cout << "point ready" << std::endl;
                    });
                }
            }
        }
        ThreadPool threads{tasks};
        threads.run(worker<typename TaskQueue::value_type>);
        
    }


    FieldValue getHarmonic(size_t ifreq, size_t irow, size_t icolumn) const 
    {
        return m_HarmonicsCatalogue[ifreq][irow][icolumn];
    }


    void constructSummary(std::vector<double> &freqs, std::vector<std::complex<double>> &amplF, const double time)
    {   
        function task{[this, &freqs, &amplF, time, w = 0, i = 0, j = 0](){}};
        using TaskQueue = ConcurentQueue<decltype(task)>;
        using ThreadPool = ThreadPoolImpl<TaskQueue>;
        
        TaskQueue tasks;
        for (size_t i = 0; i < mn_points; i++)
        {
            for (size_t j = 0; j < mn_points; j++)
            {
                for (size_t w = 0; w < mn_harm; w++)
                {   
                    tasks.push([this, &freqs, &amplF, time, w, i, j](){
                        m_SummField[i][j].E=m_SummField[i][j].E+m_HarmonicsCatalogue[w][i][j].E*amplF[w]*exp(std::complex<double>(0.0, 1.0) * freqs[w] * time);
                        //std::cout << "point ready" << std::endl;
                    });
                }
            }
        }
        ThreadPool threads{tasks};
        threads.run(worker<typename TaskQueue::value_type>);
    }


    FieldValue getSummary(size_t irow, size_t icolumn) const 
    {
        return m_SummField[irow][icolumn];
    }
private:
    ISurface& m_surf;
    SurfaceRegion m_region;
    std::vector<StrattonChuReflection> &m_fields;    
    size_t mn_harm;
    size_t mn_points;
    //std::array<std::vector<std::vector<FieldValue>>, harmonics_number> m_HarmonicsCatalogue (mn_points, std::vector<FieldValue>(mn_points));

    //std::vector<std::vector<FieldValue> m_HarmonicsCatalogue (n_points, std::vector<FieldValue>(n_points))[harmonics_number];
    
    std::vector<std::array<std::array<FieldValue, n_points>, n_points>> m_HarmonicsCatalogue;

    std::vector<std::array<FieldValue, n_points>> m_SummField;
};



#endif // FIELDSINREGION_HPP