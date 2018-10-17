#ifndef DCELL_H
#define DCELL_H

#include "gheader.h"
#include <Eigen/Dense>

class Dcell{
    typedef Eigen::Matrix3d                         cel_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,3>  pos_type;

public:
    //cell
    cel_type CLP;   // cell parameter
    
    //position
    vector<string>  Aname;  // atom name
    vector<int>     Anum;   // atom numbers
    vector<double>  Amass;  // atom mass
    vector<double>  Acharge;// atom charge
    pos_type POS;   // atom positions
    pos_type fPOS;  // fractional atom positions
    
    //wannier
    pos_type WAN;   // wannier centers
    pos_type fWAN;  // fractional wannier centers
    

    
    
    
};

#endif
