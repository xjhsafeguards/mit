#ifndef CELL_H
#define CELL_H

#include "Box.h"
#include "Position.h"

#include <vector> // std::vector
#include <string> // std::string

namespace Cell_space {
    
    //  typically box_parameters stored in Angstrom and positions in fractional
    //  periodical boundary condition
    //  have interface with
    //      position
    //      Eigen::Matrix<double,1,3,Eigen::RowMajor>
    //      vector<double>(3) // not support yet
    //      all treated as fractional
    struct Cell{
      
        typedef Eigen::Matrix<double,1,3,Eigen::RowMajor> pos_type;
        
        Box cellp;
        Positions atoms,wanniers;

        // Cell information
        int snapshot;
        double time;
        
        // atom data
        int na = -1; // num of total atoms
        std::vector<int> types_index; //type along the pos vector
        std::vector<std::string> types;
        std::vector<std::string> type_list; //list of types
        std::vector<int> number_list;
        std::vector<double> mass_list;
        
    //wrappers for Box functions
        double volume() const{return cellp.volume();}
        //  operations on positions
        const pos_type to_cart(const pos_type& ip) const{return cellp.to_cart(ip);}
        const pos_type to_frac(const pos_type& ip) const{return cellp.to_frac(ip);}
        double distance(const pos_type& ip1,const pos_type& ip2) const{return cellp.distance(ip1,ip2);}
        double angle(const pos_type& ip1,const pos_type& ip2,const pos_type& ip3) const{return cellp.angle(ip1,ip2,ip3);}
        double fdistance(const pos_type& ip1,const pos_type& ip2) const{return cellp.fdistance(ip1,ip2);}
        double fangle(const pos_type& ip1,const pos_type& ip2,const pos_type& ip3) const{return cellp.fangle(ip1,ip2,ip3);}
        // operations on atoms
        double atom_distance(int i,int j) const{return cellp.distance(atom(i),atom(j));}
        double adistance(int i,int j) const{return cellp.distance(atom(i),atom(j));}
        double atom_angle(int i,int j,int k) const{return cellp.angle(atom(i),atom(j),atom(k));}
        double aangle(int i,int j,int k) const{return cellp.angle(atom(i),atom(j),atom(k));}
        
        
    //atoms functions
        
        
    //wrappers for atoms functions
        template<typename... Ts>
        auto atom(Ts&&... ts) const -> const decltype(atoms.row(std::forward<Ts>(ts)...)){
            return atoms.row(std::forward<Ts>(ts)...);
        }
    
    //wrappers for wanniers functions
        template<typename... Ts>
        auto wannier(Ts&&... ts) const -> const decltype(wanniers.row(std::forward<Ts>(ts)...)){
            return wanniers.row(std::forward<Ts>(ts)...);
        }
    };
    
    struct Cell_qecp: public Cell{
        void read_cellp(std::istream&);
        void skip_cellp(std::istream&) const;
        void read_atoms(std::istream&);
        void skip_atoms(std::istream&) const;
        void read(std::istream&,std::istream&);
        void skip(std::istream&,std::istream&) const ;
    };
}
using namespace Cell_space;
#endif
