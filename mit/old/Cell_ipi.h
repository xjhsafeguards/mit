#ifndef CELL_IPI_H
#define CELL_IPI_H

#include "Cell.h"
#include "Physics_unit.h"

class Cell_ipi_c : public Cell{
    int cell_type;
public:
    Cell_ipi_c(int in_type=1): cell_type(in_type){is_frac=false;}
    
    virtual void read(std::istream&);
    virtual void skip(std::istream&) const ;
    
    virtual double distance(int i,int j) const {return cal_distance_c(atom_position(i),atom_position(j),cell_parameters);}
    virtual double angle(int i,int j, int k) const {return cal_angle_c(atom_position(i),atom_position(j),atom_position(k),cell_parameters);}
    virtual std::string type(int i) const {return types[i];}
    
protected:
    void read_type1(std::istream&);
    void skip_type1(std::istream&) const;
};

#endif
