#ifndef POSITION_H
#define POSITION_H

#include <Eigen/Dense>
#include <vector> // std::vector
#include <iostream> // std::istream

namespace Cell_space{
    //  typically a set of positions stored in fractional
    //  have interface with
    //      Eigen::Matrix<double,1,3,Eigen::RowMajor>
    //      vector<double>(3) // not support yet
    //      all treated as fractional
    class Positions: public Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor>{
    public:
        typedef Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> data_type;
        typedef Eigen::Matrix3d box_type;
        typedef Eigen::Matrix<double,1,3,Eigen::RowMajor> pos_type;
        typedef std::vector<double> vpos_type;
        
        //set values
        // set length of positions
        void init(int i){resize(i,3);}
        void set_n(int i){resize(i,3);}
        void load(int i,double x,double y,double z){row(i) << x,y,z;}
        void load1(int i,double x,double y,double z){row(i) << x,y,z;}
        void load1(int i,std::istream& is){double x,y,z; is>>x>>y>>z; load1(i,x,y,z);}
        // modify data
        void set_unit(double in_unit){*this *= in_unit;}
        void frac(const box_type& inbox){ *this *= inbox.inverse();}
        void set_frac(const box_type& inbox){ *this *= inbox.inverse();}
        void cart(const box_type& inbox){ *this *= inbox;}
        void set_cart(const box_type& inbox){ *this *= inbox;}
        
        //achieving position
        // const aliasing for row
        template<typename... Ts>
        auto pos(Ts&&... ts) const -> const decltype(row(std::forward<Ts>(ts)...)){
            return row(std::forward<Ts>(ts)...);
        }
        // vector<double>(3)
        const vpos_type vpos(int i) const {return to_vpos(pos(i));}
        
    private:
        const vpos_type to_vpos(const data_type& d) const{return vpos_type(d.data(),d.data()+3);}
    };
    
}

using namespace Cell_space;
#endif
