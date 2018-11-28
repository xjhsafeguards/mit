#ifndef REORGANIZE_H
#define REORGANIZE_H

#include "cellfile.h"
#include "waterfile.h"

class Reorganize
{
public:
    
    //Reorganize(){};
    //~Reorganize(){};
    
    static void Print_xyz(ofstream &ofs, Cellfile &cel);
    static void Print_xyz_water(ofstream &ofs, Cellfile &cel);
    static void Print_xyz_non_water(ofstream &ofs, Cellfile &cel);
    static void Read_average_cel(int type=0); // read average cell to INPUT.Celldm; and count to INPUT.ss_n;
    
    static void Routine();
    
};

#endif
