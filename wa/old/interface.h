#ifndef INTERFACE_H
#define INTERFACE_H

#include "density.h"
class Density;

class Interface
{
public:
    
    Interface();
    ~Interface();
    
    int snapshot;
    double time;
    int num_z;
    double *posz; // distance of atoms to the interface;
    bool allocate_posz;
    
    double get_posz(int type, int num)
    {
        assert(allocate_posz);
        return posz[type_mark[type]+num];
    }
    
    void organize_posz(const double &displacement,const double &BC)
    //organize the value of posz to a required displacement
    {
        assert(allocate_posz);
        for(int i=0; i< num_z; i++)
        {
            if(posz[i] < displacement)
                posz[i] += BC;
        }
    }
    
    /********* ILI surface ********/
    int Read_pos_ILI(ifstream &ifs, bool skip=0); //read distance to the ILI of each atom from pos_ILI file
    
    /********* Gibbs surface ********/
    double Gs; // record z value for Gs
    void Read_Gibbs_surface(const Density &dens);

private:
    int *type_mark; // mart the fisrt index of each type of atoms

};

#endif
