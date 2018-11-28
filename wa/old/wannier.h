#ifndef WANNIER_H
#define WANNIER_H

#include "water.h"

class Wannier
{
public:
    
    Wannier(){};
    ~Wannier(){};
    
    int numWC;
    vector<int> index_WC;
    Vector3<double> dipole;//(Debye) 1Debye=0.20819434eA
    // 1a.u. time =2.418884326505×10-17 s

    
};

#endif
