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
#include <memory> // shared_ptr make_shared

class Atom;
class Atom_group;
class Cell;
class Cell_manip;
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
    virtual Atom_group atoms() const;
    virtual Atom_group atoms(std::string in_type) const;
    virtual Atom_group atoms(int start_index,int end_index) const;
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
    //pos_type atom_position(int i) const {return positions.row(i);}
    pos_type atom_position(int i) const {return std::move(positions.row(i));}
    
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
    pos_type cal_frac(const pos_type& cp,const box_type& bc) const {return std::move(cp*bc.inverse());}
    pos_type cal_cart(const pos_type& fp,const  box_type& bc) const {return std::move(fp*bc);}
    
    
    // basic functions to calculate relationship between positions
    double cal_distance_f(const pos_type& fp1,const pos_type& fp2,const box_type& bc) const
    {
        return std::move(cal_cart(shortest_fvector(fp1,fp2),bc).norm());
    }
    double cal_angle_f(const pos_type& fp1,const pos_type& fp2,const pos_type& fp3,const box_type& bc) const
    {
        pos_type v1=std::move(cal_cart(shortest_fvector(fp2,fp1),bc));
        pos_type v2=std::move(cal_cart(shortest_fvector(fp3,fp1),bc));
        return std::move(acos(v1.dot(v2)/v1.norm()/v2.norm())*180/PI);
    }
    pos_type shortest_fvector(const pos_type& fp1,const pos_type& fp2) const;
    
    // basic fucntions to calculate relationship between positions for cartesian only used if box is diagonal
    double cal_distance_c(const pos_type& cp1,const pos_type& cp2,const box_type& bc) const
    {
        return std::move(shortest_vector(cp1,cp2,bc).norm());
    }
    double cal_angle_c(const pos_type& cp1,const pos_type& cp2,const pos_type& cp3,const box_type& bc) const
    {
        pos_type v1=std::move(shortest_vector(cp2,cp1,bc));
        pos_type v2=std::move(shortest_vector(cp3,cp1,bc));
        return std::move(acos(v1.dot(v2)/v1.norm()/v2.norm())*180/PI);
    }
    pos_type shortest_vector(const pos_type& cp1,const pos_type& cp2,const box_type& bc) const;
};

class Cell_manip{
    const Cell& cel;
public:
    //Cell_manip():cel(*(new Cell())){}
    Cell_manip(const Cell& in_cel):cel(in_cel){}
};

class Atom{
    //friend class Atom_group;
    const Cell& cel;
    int index;
public:
    
    typedef typename Cell::pos_type pos_type;
    
    Atom(const Cell& in_cel,int in_index):cel(in_cel),index(in_index){};
    Atom(const Atom&) = default;
    Atom(Atom&&) = default;
    
    //pos_type position() const;
    double distance(int j) const;
    double distance(const Atom& a1) const;
    double angle(int j, int k) const;
    double angle(const Atom& a1, const Atom& a2) const;
    std::string type() const;
    int get_index() const {return index;}
};

class Atom_group{
public:
    typedef std::shared_ptr<Atom> atom_type;
    typedef std::vector<std::shared_ptr<Atom> > data_type;
protected:
    const Cell& cel;
    data_type atoms;
public:
    Atom_group(const Cell& in_cel,std::vector<int> in_indexs=std::vector<int>()):cel(in_cel){
        for( const auto& index: in_indexs){
            atoms.push_back(new_atom(index));
        }
    }
    Atom_group(const Atom_group&) = default;
    Atom_group(Atom_group&&) = default;
    
    void add(const Atom& a){  atoms.push_back(new_atom(a));}
    void add(Atom&& a){  atoms.push_back(new_atom(a));}
    data_type::iterator begin() {return atoms.begin();}
    data_type::iterator end() {return atoms.end();}
    
    int size() {return atoms.size();}
protected:
    atom_type new_atom(int index) {return std::move(std::make_shared<Atom>(cel,index));}
    atom_type new_atom(const Atom& a) {return std::move(std::make_shared<Atom>(a));}
    atom_type new_atom(Atom&& a) {return std::move(std::make_shared<Atom>(a));}
};
/*
class Atom_group{
    //friend class Atom_iterator;
protected:
    const Cell& cel;
    std::vector<int> indexs;
public:
    Atom_group(const Cell& in_cel,std::vector<int> in_indexs=std::vector<int>()):cel(in_cel),indexs(in_indexs){}
    Atom_group(const Atom_group&) = default;
    Atom_group(Atom_group&&) = default;
    
    class Atom_iterator{
        const Atom_group& AG;
        std::vector<int>::iterator it;
    public:
        typedef Atom value_type;
        
        Atom_iterator(const Atom_group& in_ag,std::vector<int>::iterator in_it): AG(in_ag),it(in_it){}
        Atom operator*() {return {AG.cel,*it};}
        Atom_iterator operator++() {++it;return *this;}
        bool operator!=(const Atom_iterator& in_AI) {return it!=in_AI.it;}
    };
    
    void add(const Atom& a){ indexs.push_back(a.index);}
    void add(Atom&& a){ indexs.push_back(a.index);}
    Atom_iterator begin() {return {*this,indexs.begin()};}
    Atom_iterator end() {return {*this,indexs.end()};}
    
    int size() {return indexs.size();}
    
};*/

class Box{
    const Cell& cel;
public:
    
    Box(const Cell& in_cel):cel(in_cel){};
};

#endif
