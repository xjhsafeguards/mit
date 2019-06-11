#ifndef CELL_LAMMPS_H
#define CELL_LAMMPS_H

#include "Cell.h"

class cell_lammps: public cell{
    //0 stand for unit in metal and only xyz cellp are given
    int cell_type = 0;
    double cell_unit= 1;
    int na=0;
    std::vector<std::string> types;
    
    std::istream& _read_timestep(std::istream&);
    std::istream& _read_n(std::istream&);
    std::istream& _read_box0(std::istream&);
    std::istream& _read_pos(std::istream&);
public:
    
    cell_lammps(std::initializer_list<std::string> li_types): types(li_types){
    }
    
    virtual std::istream& read(std::istream&);
};

#endif
