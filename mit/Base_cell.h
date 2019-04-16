#ifndef BASE_CELL_H
#define BASE_CELL_H

#include <Eigen/Dense> //Matrix
#include <iostream> // istream ostream
#include <cassert> // assert
#include "Math_const.h" // PI
#include <cmath> // acos abs
#include <vector> // std::vector
#include <memory> // std::shared_ptr, std::make_shared()

namespace Cell {
    
    typedef Eigen::Matrix<double,1,3,Eigen::RowMajor> Vector;
    
    class Box;
    class Positions;
    class Position;
    class Cell;
    
    class Box : public Eigen::Matrix3d{
        typedef Eigen::Matrix3d data_type;
        typedef Eigen::Matrix3d super;
    public:
        //constructors
        using super::super;
        //IO
        inline std::istream& read(std::istream&);
        inline void read(double,double,double);
        inline std::ostream& write(std::ostream&) const;
        //Operation
        inline double volume() const;
        inline Position to_fposition(const Position& cp) const;
        inline Position to_cposition(const Position& fp) const;
        inline Vector shortest_fvector(const Position& fp1,const Position& fp2) const;
        inline Vector shortest_cvector(const Position& cp1,const Position& cp2) const;
        inline double fdistance(const Position& fp1,const Position& fp2) const;
        inline double cdistance(const Position& cp1,const Position& cp2) const;
        inline double fangle(const Position& fp1,const Position& fp2,const Position& fp3) const;
        inline double cangle(const Position& cp1,const Position& cp2,const Position& cp3) const;
    };
    
    class Positions : public Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor>{
        typedef Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> data_type;
        typedef Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> super;
    public:
        //constructors
        using super::super;
        //IO
        inline void init(int line);
        inline void read(int,double,double,double);
        inline std::istream& read(std::istream&);
        inline std::istream& read(int,std::istream&);
        inline std::ostream& write(std::ostream&) const;
        //Operations
    };
    
    class Position : public Eigen::Matrix<double,1,3,Eigen::RowMajor>{
        typedef Eigen::Matrix<double,1,3,Eigen::RowMajor> data_type;
        typedef Eigen::Matrix<double,1,3,Eigen::RowMajor> super;
    public:
        //constructors
        using super::super;
        //IO
        inline void read(double,double,double);
        inline std::istream& read(std::istream&);
        inline std::ostream& write(std::ostream&) const;
    };
    
    template <typename T>
    class Pointer_vector : public
    
    class Position_pv : public std::vector<std::shared_ptr<Position>>{
        
    };
    
    
    class Cell{
        
    };
    
    /* Box Functions */
    inline std::istream& Box::read(std::istream& is){
        assert(is.good());
        double item = 0.0;
        for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++)
            {
                is >> item;
                (*this)(row, col) = item;
            }
        return is;
    }
    inline void Box::read(double a,double b,double c){
        //this->block(0,0,3,3) = Eigen::Vector3d(a,b,c).asDiagonal();
        *this = Eigen::Vector3d(a,b,c).asDiagonal();
    }
    inline std::ostream& Box::write(std::ostream& os) const{
        assert(os.good());
        os << *this;
        return os;
    }
    inline double Box::volume() const{
        return abs(row(0).cross(row(1)).dot(row(2)));
    }
    inline Position Box::to_fposition(const Position& cp) const{
        return cp*inverse();
    }
    inline Position Box::to_cposition(const Position& fp) const{
        return fp*(*this);
    }
    inline Vector Box::shortest_fvector(const Position& fp1,const Position& fp2) const{
        Vector result = fp1 - fp2;
        auto iter = result.data();
        for(int i=0;i<3;i++){
            while(*iter<=-0.5) *iter+=1;
            while(*iter>0.5) *iter-=1;
            ++iter;
        }
        return result;
    }
    inline Vector Box::shortest_cvector(const Position& cp1,const Position& cp2) const{
        //assert()
        Vector result = cp1 - cp2;
        auto iter=result.data();
        for(int i=0;i<3;i++){
            double bcv=(*this)(i,i);
            double hbcv=bcv/2;
            while(*iter<=-hbcv) *iter+=bcv;
            while(*iter>hbcv) *iter-=bcv;
            ++iter;
        }
        return result;
    }
    inline double Box::fdistance(const Position& fp1,const Position& fp2) const{
        return to_cposition(shortest_fvector(fp1,fp2)).norm();
    }
    inline double Box::cdistance(const Position& cp1,const Position& cp2) const{
        return shortest_cvector(cp1,cp2).norm();
    }
    inline double Box::fangle(const Position& fp1,const Position& fp2,const Position& fp3) const{
        Vector v1=to_cposition(shortest_fvector(fp1,fp2));
        Vector v2=to_cposition(shortest_fvector(fp1,fp3));
        return acos(v1.dot(v2)/v1.norm()/v2.norm())*180/PI;
    }
    inline double Box::cangle(const Position& cp1,const Position& cp2,const Position& cp3) const{
        Vector v1=shortest_cvector(cp1,cp2);
        Vector v2=shortest_cvector(cp1,cp3);
        return acos(v1.dot(v2)/v1.norm()/v2.norm())*180/PI;
    }
    
    /* Positions Functions */
    inline void Positions::init(int line){
        resize(line,3);
    }
    inline void Positions::read(int i,double a,double b,double c){
        row(i) << a,b,c;
    }
    inline std::istream& Positions::read(std::istream& is){
        assert(is.good());
        double a,b,c;
        for(int __row = 0; __row < rows(); ++__row){
            is >> a >> b >> c;
            row(__row) << a,b,c;
        }
        return is;
    }
    inline std::istream& Positions::read(int i,std::istream& is){
        assert(is.good());
        double a,b,c;
        is >> a >> b >> c;
        row(i) << a,b,c;
        return is;
    }
    inline std::ostream& Positions::write(std::ostream& os) const{
        os << *this;
        return os;
    }
    
    /*
     class Position: public Eigen::Matrix<double,1,3,Eigen::RowMajor>{
     typedef Eigen::Matrix<double,1,3,Eigen::RowMajor> data_type;
     typedef Eigen::Matrix<double,1,3,Eigen::RowMajor> super;
     public:
     //constructors
     using super::super;
     //IO
     inline std::istream& read(std::istream&);
     inline void read(double,double,double);
     inline std::ostream& write(std::ostream&) const;
     //Operations
     Position to_fraction() const;
     
     };*/

}

using namespace Cell;
#endif
