#ifndef DCELL_H
#define DCELL_H

#include "gheader.h"
#include <Eigen/Dense>

class Dcell{
    typedef Eigen::Matrix3d                         cel_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,3>  nvec_type;
    typedef Eigen::RowVector3d                      vec_type;

public:
    //cell
    cel_type CLP;   // cell parameter
    
    //position
    vector<int>     Anum;   // atom numbers
    vector<string>  Aname;  // atom name
    vector<double>  Amass;  // atom mass
    vector<double>  Acharge;// atom charge
    nvec_type POS;   // atom positions
    nvec_type fPOS;  // fractional atom positions
    
    //wannier
    nvec_type WAN;   // wannier centers
    nvec_type fWAN;  // fractional wannier centers
    
    //I/O
    void set_cel(cel_type &CLP_in);
    void set_pos(nvec_type &POS_in,vector<int> Anum_in,vector<string> Aname_in=vector<string>(),vector<double> Amass_in=vector<double>(),vector<double> Acharge_in=vector<double>());
    

    
    
    
};

#endif
