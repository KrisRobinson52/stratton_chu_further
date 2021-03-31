#ifndef CSVSAVER_HPP
#define CSVSAVER_HPP

#include "stratton-chu/types.hpp"
#include "stratton-chu/field.hpp"

#include <fstream>

#include <string>

class CSVSaver
{
public:
    CSVSaver(const char* filename);
    void add_point(Position p, FieldValue fv);

private:
    std::ofstream m_ofs;
};

#endif // CSVSAVER_HPP
