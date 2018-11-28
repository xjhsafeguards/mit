#ifndef WATERFILE_H
#define WATERFILE_H

#include "cellfile.h"
#include "water.h"

class Waterfile
{
public:
    
    Waterfile();
    ~Waterfile();
    
    int O,H,C; //mark the index of O and H
    
    int nwater; //total number of water(include ion)
    Water *waters; //holds the information of water
    
    int nion;  //total number of ions
    vector<int> ions;  //mark the index of ions
    
    bool allocate_waters;
    const Cellfile *cel;
    
    void Read_water(Cellfile &Cel); //get water from given cellfile
    //void Mark_water(Cellfile &Cel,string tmpid = "S"); // use different id to mark water
    
    
    static void Print_water(string out_file = "Water.txt");
    static void Add_constrain();
    
    
public:
    
    Vector3<double> posO(int nO) const; //return the position of nth O in Cellfile unit
    Vector3<double> posH(int nO, int nH) const; //return the position of indexth water's indexth H in Cellfile unit
    static void Animate_water();
    
    //add 2018.3.30
    int Show_nwater() const {return nwater;}
    int Show_nH(int nO) const {return waters[nO].Show_nH();}
    int Show_indexH(int nO, int nH) const {return waters[nO].Show_indexH(nH);}
    
};

#endif
