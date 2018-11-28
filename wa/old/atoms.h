#ifndef ATOMS_H
#define ATOMS_H

#include "vec3.h"
#include "gfun.h"


class Atoms
{
    public:
    
    Atoms();
    ~Atoms();
    
    string id;  //atom id
    int na; //total num of this Atoms
    Vector3<double> *pos; //the postion of atoms(input)
    
    int *index; //the index in the in_file
    
    int snapshot;
    
    double mass;//(SI)
    double charge;//(SI)
    
    bool allocate_pos;
    bool allocate_index;
    
    public:
    
    void read_pos(ifstream &ifs, bool frac=false); //read the fractional coordinates or cartesian coordinates from input file
    void read_pos2(ifstream &ifs, bool frac=false);// read given type of atoms from input file
    
    void read_geo(ifstream &ifs);//read cartesian position from .pos file within one snapshot
    
};

#endif
