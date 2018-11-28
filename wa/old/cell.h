#ifndef CELL_H
#define CELL_H

#include "atoms.h"

class Cell
{
public:
    
    Cell(){};
    ~Cell(){};
    
    int snapshot;
    double time;
    
    int ntype;
    int nband;
    int natom;
    Vector3<double> celldm;//(input)
    
    //information for each atom
    Atoms *atoms;// saving the atom information in a cell
    Vector3<double> *wan_centers;//(input)
};

#endif
