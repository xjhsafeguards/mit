#ifndef CELL_H
#define CELL_H

#include <Eigen/Dense>
#include "Math_const.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

class Atom;
class Atom_group;
class Cell;
class Box;


class Cell{
public:
    
    friend class Atom;
    friend class Box;
    typedef Eigen::Matrix3d box_type;
    typedef Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> posf_type;
    typedef Eigen::Matrix<double,1,3,Eigen::RowMajor> pos_type;
    
    Cell(){}
    
    bool operator==(const Cell& cel) const {return this==&cel;}
    
    virtual Atom atom(int index) const;
    virtual Box  box() const;
    
    // Reading in Cells
    virtual void read(std::istream&){std::cerr << "no read implemented! Quit!"; std::exit(0);}
    virtual void skip(std::istream&) const {std::cerr << "no skip implemented! Quit!"; std::exit(0);}
    virtual void read(std::istream&,std::istream&){std::cerr << "no read implemented! Quit!"; std::exit(0);}
    virtual void skip(std::istream&,std::istream&) const {std::cerr << "no skip implemented! Quit!"; std::exit(0);}
    virtual void read_cell(std::istream&){std::cerr << "no read cell implemented! Quit!"; std::exit(0);}
    virtual void skip_cell(std::istream&) const {std::cerr << "no skip cell implemented! Quit!"; std::exit(0);}
    virtual void read_pos(std::istream&){std::cerr << "no read pos implemented! Quit!"; std::exit(0);}
    virtual void skip_pos(std::istream&) const {std::cerr << "no read pos implemented! Quit!"; std::exit(0);}
    // set Cells
    virtual void set_box(double a, double b, double c) {read_box(a,b,c);}
    
    // Writing Cells
    virtual std::ostream& write_cell(std::ostream& os) const { os << cell_parameters; return os;}
    virtual std::ostream& write_positions(std::ostream& os) const { os << positions; return os;}
    virtual std::ostream& write_position(std::ostream& os, int i) const { os << atom_position(i); return os;}
    
    virtual double volume() const;
    virtual double snapshot() const {return c_snapshot;}
    virtual double time() const {return c_time;}
    virtual std::string type(int) const {std::cerr << "no type implemented! Quit!"; std::exit(0);}
    //distance between atom i and atom j;
    virtual double distance(int i,int j) const;
    //angle of atom j-i-k
    virtual double angle(int i,int j, int k) const;
    
    //density of the system
    //virtual double density(){std::cerr << "no density implemented! Quit!"; std::exit(0);}
    
protected:
    
    // data
    bool is_frac=true;
    int c_snapshot=-1;
    double c_time=-1;
    box_type cell_parameters;
    posf_type positions;
    
    // atom data
    int na; // num of total atoms
    std::vector<int> types_index; //type along the pos vector
    std::vector<std::string> types;
    std::vector<std::string> type_list; //list of types
    std::vector<int> number_list;
    std::vector<double> mass_list;
    
    //get atom_index
    //int atom_index(std::string name, int index){}
    
    //get position
    pos_type atom_position(int i) const {return positions.row(i);}

//read
    // init positions for read
    void init_positions(int i) {positions.resize(i,3);}
    void read_positions(int i, double a, double b, double c) {positions.row(i) << a,b,c;}
    void read_box(double a, double b, double c) {cell_parameters = Eigen::Vector3d(a,b,c).asDiagonal();}
    
//unit
    void set_positions_unit(double unit){positions*=unit;}
    void set_box_unit(double unit){cell_parameters*=unit;}
    void set_unit(double unit){positions*=unit;cell_parameters*=unit;}
    
    // convert cell to fractional or cartesian
    void to_frac();
    void to_cart();
    // basic functions to modify positions
    pos_type cal_frac(const pos_type& cp,const box_type& bc) const {return cp*bc.inverse();}
    pos_type cal_cart(const pos_type& fp,const  box_type& bc) const {return fp*bc;}
    
    
    // basic functions to calculate relationship between positions
    double cal_distance_f(const pos_type& fp1,const pos_type& fp2,const box_type& bc) const
    {
        return cal_cart(shortest_fvector(fp1,fp2),bc).norm();
    }
    double cal_angle_f(const pos_type& fp1,const pos_type& fp2,const pos_type& fp3,const box_type& bc) const
    {
        pos_type v1=cal_cart(shortest_fvector(fp2,fp1),bc);
        pos_type v2=cal_cart(shortest_fvector(fp3,fp1),bc);
        return acos(v1.dot(v2)/v1.norm()/v2.norm())*180/PI;
    }
    pos_type shortest_fvector(const pos_type& fp1,const pos_type& fp2) const;
    
    // basic fucntions to calculate relationship between positions for cartesian only used if box is diagonal
    double cal_distance_c(const pos_type& cp1,const pos_type& cp2,const box_type& bc) const
    {
        return shortest_vector(cp1,cp2,bc).norm();
    }
    double cal_angle_c(const pos_type& cp1,const pos_type& cp2,const pos_type& cp3,const box_type& bc) const
    {
        pos_type v1=shortest_vector(cp2,cp1,bc);
        pos_type v2=shortest_vector(cp3,cp1,bc);
        return acos(v1.dot(v2)/v1.norm()/v2.norm())*180/PI;
    }
    pos_type shortest_vector(const pos_type& cp1,const pos_type& cp2,const box_type& bc) const;
};


class Atom{
    friend class Atom_group;
    const Cell& cel;
    int index;
public:
    
    typedef typename Cell::pos_type pos_type;
    
    Atom(const Cell& in_cel,int in_index):cel(in_cel),index(in_index){};
    
    //pos_type position() const;
    double distance(int j) const;
    double distance(const Atom& a1) const;
    double angle(int j, int k) const;
    double angle(const Atom& a1, const Atom& a2) const;
    std::string type() const;
};

class Atom_group{
    const Cell& cel;
    std::vector<int> indexs;
public:
    Atom_group(const Cell& in_cel,std::vector<int> in_indexs=std::vector<int>()):cel(in_cel),indexs(in_indexs){}
    
    void add(const Atom& a){ indexs.push_back(a.index);}
};

class Box{
    const Cell& cel;
public:
    
    Box(const Cell& in_cel):cel(in_cel){};
};

#endif
