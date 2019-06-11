#ifndef BOX_H
#define BOX_H

#include <Eigen/Dense> //Eigen::Matrix3d
#include <iostream> //std::istream
#include <cmath> //std::abs std::remainder
#include <vector> //std::vector

#include "Math_const.h"

namespace Cell_space{
    //  typically box_parameters stored in Angstrom
    //  periodical boundary condition
    //  have interface with
    //      position
    //      Eigen::Matrix<double,1,3,Eigen::RowMajor>
    //      vector<double>(3) // not support yet
    //      all treated as fractional
    class Box: public Eigen::Matrix3d{
    public:
        typedef Eigen::Matrix3d data_type;
        typedef Eigen::Matrix<double,1,3,Eigen::RowMajor> pos_type;
        typedef Eigen::Matrix3d super;
        
        //constructors
        using super::super;
        
        //set values
        void load(std::istream&);
        void load(double, double, double);
        void load(double, double, double, double, double, double, double, double, double);
        void load3(std::istream&);
        void load3(double, double, double);
        //modify values
        void set_unit(double);
        
        //operations
        double volume() const;
        
        //operations on positions
        const pos_type to_cart(const pos_type&) const;
        const pos_type to_frac(const pos_type&) const;
        double distance(const pos_type&,const pos_type&) const;
        double angle(const pos_type&,const pos_type&,const pos_type&) const;
        double fdistance(const pos_type&,const pos_type&) const;
        double fangle(const pos_type&,const pos_type&,const pos_type&) const;
        
        //interface to vector<double>(3)
        double distance(const std::vector<double>&,const std::vector<double>&);
        double angle(const std::vector<double>&,const std::vector<double>&,const std::vector<double>&);
        
    private:
        //returns the shortest vector considering boundary condition;
        const pos_type _shortest_vector(const pos_type&,const pos_type&) const;
    };
    
    inline void Box::load(std::istream& is){
        double a[9];
        auto it=a;
        for(int i=0; i<9; ++i){
            is >> *it;
            ++it;
        }
        *this=Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor>>(a);
    }
    inline void Box::load(double a, double b, double c){
        load3(a,b,c);
    }
    inline void Box::load(double a, double b, double c, double d, double e, double f, double g, double h, double i){
        *this << a,b,c,d,e,f,g,h,i;
    }
    inline void Box::load3(std::istream& is){
        double a,b,c;
        is >> a >> b >> c;
        *this << a,0,0,0,b,0,0,0,c;
    }
    inline void Box::load3(double a, double b, double c){
        *this << a,0,0,0,b,0,0,0,c;
    }
    inline void Box::set_unit(double in_unit){
        *this *= in_unit;
    }
    inline double Box::volume() const{
        return std::abs(row(0).cross(row(1)).dot(row(2)));
    }
    inline const Box::pos_type Box::to_cart(const pos_type& p) const{
        return p*(*this);
    }
    inline const Box::pos_type Box::to_frac(const pos_type& p) const{
        return p*(this->inverse());
    }
    inline double Box::distance(const pos_type& p1,const pos_type& p2) const{
        return fdistance(p1,p2);
    }
    inline double Box::angle(const pos_type& p1,const pos_type& p2,const pos_type& p3) const{
        return fangle(p1,p2,p3);
    }
    inline double Box::fdistance(const pos_type& p1,const pos_type& p2) const{
        return to_cart(_shortest_vector(p1,p2)).norm();
    }
    inline double Box::fangle(const pos_type& p1,const pos_type& p2,const pos_type& p3) const{
        pos_type v1=to_cart(_shortest_vector(p2,p1));
        pos_type v2=to_cart(_shortest_vector(p3,p1));
        return std::move(acos(v1.dot(v2)/v1.norm()/v2.norm())*180/PI);
    }
    
    inline const Box::pos_type Box::_shortest_vector(const pos_type& p1,const pos_type& p2) const{
        return pos_type(std::remainder(p2[0]-p1[0],1),std::remainder(p2[1]-p1[1],1),std::remainder(p2[2]-p1[2],1));
    }
}

using namespace Cell_space;
#endif
