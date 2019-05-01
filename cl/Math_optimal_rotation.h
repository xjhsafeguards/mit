#ifndef MATH_OPTIMAL_ROTATION
#define MATH_OPTIMAL_ROTATION

#include <vector>
#include <Eigen/Dense>
#include "Math_linearalgebra.h"

//********************************************************//
// A single function wrapper for solving rotation problem
// Read in two n*3 vectors
// Output the rotation matrix from v1 to v2
//********************************************************//

class Optimal_rotation{
    bool _allow_reflection=false;
    Eigen::MatrixXd M1,M2,MR;
    Eigen::Matrix3d H,R;
    void compute_covariance_matrix();
    void compute_rotation_matrix();
    void compute_rotated_vectors();
    
public:
    Optimal_rotation(){};
    Optimal_rotation(const std::vector <std::vector<double> >& v1,const std::vector <std::vector<double> >& v2);
    ~Optimal_rotation(){};
    
    void allow_reflection(){_allow_reflection=true;}
    void no_reflection(){_allow_reflection=false;}
    std::vector <std::vector<double> > Solve(const std::vector <std::vector<double> >& v1,const std::vector <std::vector<double> >& v2);
    std::vector <std::vector<double> > Solve(const std::vector <std::vector<double> >& v1,const std::vector <std::vector<double> >& v2, std::vector <std::vector<double> >& Rotation_mat);
    std::vector <Vector3<double> > Solve(const std::vector <Vector3<double> >& v1,const std::vector <Vector3<double> >& v2);
    std::vector <std::vector<double> > Rotation_matrix();
    std::vector <std::vector<double> > Rotated_vectors();
};

/*
class Optimal_rotation{
    void compute_covariance_matrix(const std::vector <std::vector<double> >& v1,const std::vector <std::vector<double> >& v2, Eigen::Matrix3d& H);
    void compute_covariance_matrix(const Eigen::MatrixXd& M1,const Eigen::MatrixXd& M2, Eigen::Matrix3d& H);
    void compute_rotation_matrix(const Eigen::Matrix3d& H, Eigen::Matrix3d& R);
public:
    Optimal_rotation(){};
    ~Optimal_rotation(){};
    void Solve(std::vector <std::vector<double> >& v1, std::vector <std::vector<double> >& v2,  std::vector <std::vector<double> >& result);
};
 */

#endif
