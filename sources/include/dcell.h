#ifndef DCELL_H
#define DCELL_H

#include <Eigen/Dense>

class Dcell{
    typedef Eigen::Matrix3d                         cel_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,3>  pos_type;

public:
    cel_type CLP;   // cell parameter
    pos_type POS;   // atom positions
    pos_type WAN;   // wannier centers
    
};

#endif
