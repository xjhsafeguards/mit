#ifndef DCELL_H
#define DCELL_H

#include "gheader.h"
#include <Eigen/Dense>

class Dcell{
    typedef Eigen::Matrix3d                         cel_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,3>  pos_type;
    typedef Eigen::RowVector3d                      vec_type;

public:
    //cell
    cel_type CLP;   // cell parameter
    
    //position
    vector<int>     Anum;   // atom numbers
    vector<string>  Aname;  // atom name
    vector<double>  Amass;  // atom mass
    vector<double>  Acharge;// atom charge
    pos_type POS;   // atom positions
    pos_type fPOS;  // fractional atom positions
    
    //wannier
    pos_type WAN;   // wannier centers
    pos_type fWAN;  // fractional wannier centers
    
    //I/O
    void set_cel(cel_type &CLP_in);
    void set_pos(pos_type &POS_in,vector<int> Anum_in,vector<string> Aname_in=vector<string>(),vector<double> Amass_in=vector<double>(),vector<double> Acharge_in=vector<double>())
    

    
    
    
};

#endif
