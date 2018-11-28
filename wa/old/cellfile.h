#ifndef CELLFILE_H
#define CELLFILE_H

#include "cell.h"

class Cellfile: public Cell
{
public:
    
    Cellfile();
    ~Cellfile();
    
    bool allocate_atoms;
    bool allocate_wan_centers;
    
    int ReadGeometry(ifstream &ifs,bool skip=0);//creat and read(or skip) data from ifs to cellfile
    void ReadInfile(ifstream &ifs);// read cellfile from the input file
    
    int ReadWC(ifstream &ifs,bool skip=0);//wannier center
    
    void organize_pos();// put pos of atoms in [0,celldm]
    void organize_wan();
    void set_celldm();  // set celldm
    
    //add 2018.4.11 by jianhang
    double Show_cell_volume() const; // in iput unit
};

#endif
