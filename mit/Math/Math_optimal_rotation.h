#ifndef MATH_OPTIMAL_ROTATION
#define MATH_OPTIMAL_ROTATION

#include <vector>
#include <Eigen/Dense>
#include "Math_linearalgebra.h"

//********************************************************//
// A single function wrapper for solving rotation problem
// Read in two n*3 vectors
//  have interface with vector<vector<double>> or vector<Vector3<double>>
// Rotation_matrix();
//  Output the rotation matrix from v1 to v2
// Rotated_vectors();
//  Output the rotated matrix of v1 to v2
//********************************************************//

class Optimal_rotation{
    bool _allow_reflection=false;
    Eigen::MatrixXd M1,M2,MR;
    Eigen::Matrix3d H,R;
    void compute_covariance_matrix();
    void compute_rotation_matrix();
    void compute_rotated_vectors();
    Eigen::MatrixXd compute_rotated_vectors(const Eigen::MatrixXd&);
    
    //format trans
    std::vector <std::vector<double> > Mtovv(const Eigen::MatrixXd&);
    std::vector <Vector3<double> > MtovV3(const Eigen::MatrixXd&);
    Eigen::MatrixXd to_M(const std::vector <std::vector<double> >&);
    Eigen::MatrixXd to_M(const std::vector <Vector3<double> >&);
    void vvtoM(const std::vector <std::vector<double> >&, Eigen::MatrixXd&);
    void vV3toM(const std::vector <Vector3<double> >&, Eigen::MatrixXd&);
    
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
    
    //different rotation use the solved Rotation_matrix
    std::vector <Vector3<double> > Rotated_vectors(const std::vector <Vector3<double> >&);
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
