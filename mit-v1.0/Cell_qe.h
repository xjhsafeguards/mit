#ifndef CELL_QE_H
#define CELL_QE_H

#include <cassert>
#include <vector>
#include <string>
#include <initializer_list>
#include <numeric>      // std::accumulate

#include "Cell.h"
#include "Physics_unit.h"

class Cell_qecp_c: public Cell{
    
public:
    Cell_qecp_c(std::initializer_list<int> li_na,std::initializer_list<std::string> li_types){
        type_list = std::vector<std::string>(li_types);
        number_list = std::vector<int>(li_na);
        assert(type_list.size()==number_list.size());
        is_frac=false;
        na = std::accumulate(number_list.cbegin(),number_list.cend(),0);
        for(int i=0; i!=type_list.size(); ++i){
            for(int j=0; j!=number_list[i]; ++j)
                types.push_back(type_list[i]);
        }
    }
    
    virtual void read_cell(std::istream&);
    virtual void skip_cell(std::istream&) const;
    virtual void read_pos(std::istream&);
    virtual void skip_pos(std::istream&) const;
    virtual void read(std::istream&,std::istream&);
    virtual void skip(std::istream&,std::istream&) const ;
    
    virtual double distance(int i,int j) const {return cal_distance_c(atom_position(i),atom_position(j),cell_parameters);}
    virtual double angle(int i,int j, int k) const {return cal_angle_c(atom_position(i),atom_position(j),atom_position(k),cell_parameters);}
    virtual std::string type(int i) const {
        return types[i];
    }
};


#endif
