#ifndef H_BONDSFILE_H
#define H_BONDSFILE_H

#include "H_bond.h"
#include "waterfile.h"

class Hbondfile: public Waterfile
{
public:
    
    Hbondfile();
    ~Hbondfile();
    
    //some basic functions
    static void Routine();
    void Read_Hbond(Cellfile &Cel);
    void Print_average(ofstream &ofs);
    void Print_each(ofstream &ofs);
    
    //save Hbonds information for each water
    Hbond *Hbonds;
    bool allocate_Hbonds;
    
    //Hbonds information for those water with numH != 2
    double avg_aion;
    double avg_dion;
    double avg_count_ion;
    
    //Hbonds information for those water with numH == 2
    double avg_awater;
    double avg_dwater;
    double avg_count_water;
    
private:
    void Donate(Cellfile &Cel,int O1, int O2);
    // record if water with index(O1) donate a H_bond to water with index (O2);
};

#endif

