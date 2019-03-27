#ifndef CELL_H
#define CELL_H

#include <Eigen/Dense>
#include "Math_const.h"
#include <cmath>
#include <iostream>
#include <stdlib.h>

class Atom;
class Cell;
class Box;

class Atom{
    
    const Cell& cel;
    int index;
    
public:
    
    Atom(const Cell& in_cel,int in_index):cel(in_cel),index(in_index){};
    double distance() const;
};

class Box{
    
    const Cell& cel;
    
public:
    
    Box(const Cell& in_cel):cel(in_cel){};
};

class Cell{
public:
    
    typedef Eigen::Matrix3d box_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> posf_type;
    typedef Eigen::Matrix<double,1,3,Eigen::RowMajor> pos_type;
    
    virtual Atom atom(int index) const {return Atom(*this,index);}
    virtual Box  box() const {return Box(*this);}
    
    // Reading in Cells
    virtual void read(std::istream&){std::cerr << "no read implemented! Quit!"; std::exit(0);}
    virtual void read(std::istream&,std::istream&){std::cerr << "no read implemented! Quit!"; std::exit(0);}
    virtual void read_cell(std::istream&){std::cerr << "no read implemented! Quit!"; std::exit(0);}
    virtual void read_pos(std::istream&){std::cerr << "no read implemented! Quit!"; std::exit(0);}
    
    //distance between atom i and atom j;
    virtual double distance() const;
    //angle of atom i-j-k
    virtual double angle() const;
    
//protected:
    
    // data
    bool is_frac=true;
    int snapshot=-1;
    double time=-1;
    box_type cell_parameters;
    posf_type positions;
    // convert cell to fractional or cartesian
    void to_frac();
    void to_cart();
    
    
    // basic functions to modify positions
    pos_type cal_frac(const pos_type& fp,const box_type& bc) const {return fp*bc.inverse();}
    pos_type cal_cart(const pos_type& cp,const  box_type& bc) const {return cp*bc;}
    // basic functions to calculate relationship between positions
    double cal_distance(const pos_type& fp1,const pos_type& fp2,const box_type& bc) const
    {
        return cal_cart(shortest_fvector(fp1,fp2),bc).norm();
    }
    double cal_angle(const pos_type& fp1,const pos_type& fp2,const pos_type& fp3,const box_type& bc) const
    {
        pos_type v1=cal_cart(shortest_fvector(fp1,fp2),bc);
        pos_type v2=cal_cart(shortest_fvector(fp3,fp2),bc);
        return acos(v1.dot(v2)/v1.norm()/v2.norm())*180/PI;
    }
    pos_type shortest_fvector(const pos_type& fp1,const pos_type& fp2) const;
    // basic fucntions to 
};

#endif
