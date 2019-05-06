#include "Math_optimal_rotation.h"
#include <Eigen/Dense>
#include <cassert>

Optimal_rotation::Optimal_rotation(const std::vector <std::vector<double> >& v1,const std::vector <std::vector<double> >& v2){
    assert(v1.size()==v2.size());
    vvtoM(v1,M1);
    vvtoM(v2,M2);
    compute_covariance_matrix();
    compute_rotation_matrix();
    compute_rotated_vectors();
}
std::vector <std::vector<double> > Optimal_rotation::Solve(const std::vector <std::vector<double> >& v1,const std::vector <std::vector<double> >& v2){
    assert(v1.size()==v2.size());
    vvtoM(v1,M1);
    vvtoM(v2,M2);
    compute_covariance_matrix();
    compute_rotation_matrix();
    compute_rotated_vectors();
    return Rotated_vectors();
}

std::vector <std::vector<double> > Optimal_rotation::Solve(const std::vector <std::vector<double> >& v1,const std::vector <std::vector<double> >& v2, std::vector <std::vector<double> >& Rotation_mat){
    assert(v1.size()==v2.size());
    vvtoM(v1,M1);
    vvtoM(v2,M2);
    compute_covariance_matrix();
    compute_rotation_matrix();
    compute_rotated_vectors();
    Rotation_mat = Rotation_matrix();
    return Rotated_vectors();
}
std::vector<Vector3<double> > Optimal_rotation::Solve(const std::vector <Vector3<double> >& v1,const std::vector <Vector3<double> >& v2){
    assert(v1.size()==v2.size());
    vV3toM(v1,M1);
    vV3toM(v2,M2);
    compute_covariance_matrix();
    compute_rotation_matrix();
    compute_rotated_vectors();
    return MtovV3(MR);
}

void Optimal_rotation::compute_covariance_matrix(){
    assert(M1.rows()==3);
    assert(M2.rows()==3);
    H = M1 * M2.transpose();
}

void Optimal_rotation::compute_rotation_matrix(){
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU| Eigen::ComputeFullV);
    R = svd.matrixV() * svd.matrixU().transpose();
    if(!_allow_reflection and R.determinant()<0){
        Eigen::Matrix3d V=svd.matrixV();
        V.col(2) *= -1;
        //R.col(2) *= -1;
        R = V*svd.matrixU().transpose();
    }
}

void Optimal_rotation::compute_rotated_vectors(){
    MR = R*M1;
}

Eigen::MatrixXd Optimal_rotation::compute_rotated_vectors(const Eigen::MatrixXd& M_in){
    return R*M_in;
}

std::vector <std::vector<double> > Optimal_rotation::Mtovv(const Eigen::MatrixXd& M_in){
    assert(M_in.rows()==3);
    std::vector <std::vector<double> > result;
    for(int i=0;i<M_in.cols();++i){
        result.push_back(std::vector<double>({M_in(0,i),M_in(1,i),M_in(2,i)}));
    }
    return result;
}

std::vector <Vector3<double> > Optimal_rotation::MtovV3(const Eigen::MatrixXd& M_in){
    assert(M_in.rows()==3);
    std::vector<Vector3<double> > result;
    for(int i=0;i<M_in.cols();++i){
        result.push_back(Vector3<double>(M_in(0,i),M_in(1,i),M_in(2,i)));
    }
    return result;
}

Eigen::MatrixXd Optimal_rotation::to_M(const std::vector <std::vector<double> >& v1){
    Eigen::MatrixXd tmp;
    vvtoM(v1,tmp);
    return tmp;
}

Eigen::MatrixXd Optimal_rotation::to_M(const std::vector <Vector3<double> >& v1){
    Eigen::MatrixXd tmp;
    vV3toM(v1,tmp);
    return tmp;
}

void Optimal_rotation::vvtoM(const std::vector <std::vector<double> >& v1, Eigen::MatrixXd& M){
    M.resize(3,v1.size());
    for(int i=0; i<v1.size(); ++i){
        M(0,i) = v1.at(i)[0];
        M(1,i) = v1.at(i)[1];
        M(2,i) = v1.at(i)[2];
    }
}
void Optimal_rotation::vV3toM(const std::vector <Vector3<double> >& v1, Eigen::MatrixXd& M){
    M.resize(3,v1.size());
    for(int i=0; i<v1.size(); ++i){
        M(0,i) = v1.at(i)[0];
        M(1,i) = v1.at(i)[1];
        M(2,i) = v1.at(i)[2];
    }
}

std::vector <std::vector<double> > Optimal_rotation::Rotation_matrix(){
    std::vector <std::vector<double> > result;
    for(int i=0;i<3;++i){
        result.push_back(std::vector<double>({R(i,0),R(i,1),R(i,2)}));
    }
    return result;
}
std::vector <std::vector<double> > Optimal_rotation::Rotated_vectors(){
    return Mtovv(MR);
}

std::vector <Vector3<double> > Optimal_rotation::Rotated_vectors(const std::vector <Vector3<double> >&v1){
    Eigen::MatrixXd _M = R*to_M(v1);
    return MtovV3(_M);
}
        

/*
void Optimal_rotation::Solve(std::vector <std::vector<double> >& v1, std::vector <std::vector<double> >& v2,  std::vector <std::vector<double> >& result){
    assert(v1.size()==v2.size());
    
    Eigen::Matrix3d H,R;
    Eigen::MatrixXd M1(3,v1.size()),M2(3,v2.size());
    
    for(int i=0; i<v1.size(); ++i){
        M1(0,i) = v1.at(i).at(0);
        M1(1,i) = v1.at(i).at(1);
        M1(2,i) = v1.at(i).at(2);
        M2(0,i) = v2.at(i).at(0);
        M2(1,i) = v2.at(i).at(1);
        M2(2,i) = v2.at(i).at(2);
    }
    compute_covariance_matrix(M1,M2,H);
    compute_rotation_matrix(H,R);
    
    result.clear();
    for(int i=0;i<3;++i){
        result.push_back(std::vector<double>({R(i,0),R(i,1),R(i,2)}));
    }
}

void Optimal_rotation::compute_covariance_matrix(const std::vector <std::vector<double> >& v1,const std::vector <std::vector<double> >& v2, Eigen::Matrix3d& H){
    H = Eigen::Matrix3d::Zero();
    auto it1=v1.cbegin();
    auto it2=v2.cbegin();
    for(;it1!=v1.cend();++it1,++it2){
        Eigen::Vector3d vec1(it1->at(0),it1->at(1),it1->at(2));
        Eigen::RowVector3d vec2(it2->at(0),it2->at(1),it2->at(2));
        H += vec1 * vec2;
    }
}

void Optimal_rotation::compute_covariance_matrix(const Eigen::MatrixXd& M1,const Eigen::MatrixXd& M2, Eigen::Matrix3d& H){
    assert(M1.rows()==3);
    assert(M2.rows()==3);
    H = M1 * M2.transpose();
}

void Optimal_rotation::compute_rotation_matrix(const Eigen::Matrix3d& H, Eigen::Matrix3d& R){
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU| Eigen::ComputeFullV);
    R = svd.matrixV() * svd.matrixU().transpose();
    if(R.determinant()<0)
        R.col(2) *= -1;
}

void Optimal_rotation::compute_rotated_vectors(const Eigen::MatrixXd& M1, const Eigen::Matrix3d& R){
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU| Eigen::ComputeFullV);
    R = svd.matrixV() * svd.matrixU().transpose();
    if(R.determinant()<0)
        R.col(2) *= -1;
}
*/
