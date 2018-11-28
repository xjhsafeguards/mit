#ifndef DENSITY_H
#define DENSITY_H

#include "cellfile.h"
#include "interface.h"
class Interface;

class Density
{
    public:
    
    Density();
    ~Density();
    
    static void Routine();
    void density_simple(ifstream &ifs, ofstream &logofs, bool print=0);
    void density_boundary(ifstream &ifs, ofstream &ofs,double displacement,bool print=0);
    
    
    int delta;
    
    bool allocate_dens;
    
    // double *mass; //to record the mass
    double *dens; //to record density per peice (SI)
    double delz;  // z length per peice (SI)
    double vol;   //volume per peice (SI)
    double *position; // to record the position per peice (SI)
    
    // print density
    void print_density(string out_file = "density.txt");
    void print_density(ofstream &ofs)
    {
        assert(ofs.good());
        for(int k=0;k<delta;k++)
        {
            ofs << position[k] << "\t" << dens[k] << endl;
        }
    }
    
    //geting density
    void allocate(double pos0=0,double pos1=0);//to allocate positions and dens according to INPUT
    void Read_cel_simple(Cellfile &cel);    //read density from cellfile to allocated densityfile
    void Read_ILI_simple(Interface &iff);   //read density from pos_ILI file
    
    
    void add(Density &tmp)
    {
        assert(delz==tmp.delz);
        for(int i=0;i<delta;i++)
        {
            dens[i] = dens[i] + tmp.dens[i];
        }
    }
    
};

#endif
