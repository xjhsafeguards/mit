#ifndef CELL_H
#define CELL_H

#include "Box.h"
#include "Position.h"
#include "Physics_unit.h"

#include <vector> // std::vector
#include <string> // std::string
#include <iostream> // std::istream
#include <sstream> // std::ostringstream

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
        int cell_type;
        int snapshot;
        double time;
        
        // atom data
        int na = -1; // num of total atoms
        std::vector<int> types_index; //type along the pos vector
        std::vector<std::string> types;
        std::vector<std::string> type_list; //list of types
        std::vector<int> number_list;
        std::vector<double> mass_list;
        
        //IO methods
        
        // read in .cel file
        void read_cel(std::istream&);
        void skip_cel(std::istream&) const;
        // read in .pos file and keep cart()
        void read_pos(std::istream&);
        void skip_pos(std::istream&) const;
        
        // read in .xyz file and returns the comment line
        std::string read_xyz(std::istream&,double in_unit=1);
        void skip_xyz(std::istream&) const;
        
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
        void atom_frac(){atoms.frac(cellp);}
        void afrac(){atoms.frac(cellp);}
        void atom_cart(){atoms.cart(cellp);}
        void acart(){atoms.cart(cellp);}
        
        
        
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
        //cellp
        void read_cellp(std::istream& is){read_cel(is);}
        void skip_cellp(std::istream& is) const{skip_cel(is);}
        //atoms
        void read_atoms(std::istream& is){read_pos(is);}
        void skip_atoms(std::istream& is) const{skip_pos(is);}
        //all
        void read(std::istream& is1,std::istream& is2){
            read_cellp(is1);
            read_atoms(is2);
            afrac();
        }
        void skip(std::istream& is1,std::istream& is2) const {
            skip_cellp(is1);
            skip_atoms(is2);
        }
    };
    
    struct Cell_ipi: public Cell{
        //cellp
        void read_cellp1(std::istream&,double in_unit=1);
        void skip_cellp1(std::istream&){}
        void read_cellp1(std::string in_string,double in_unit=1){std::istringstream is_com(in_string);read_cellp1(is_com,in_unit);}
        void skip_cellp1(std::string){}
        //atoms
        void read_atoms(std::istream& is,double in_unit=1){read_xyz(is,in_unit);}
        void skip_atoms(std::istream& is) const{skip_xyz(is);}
        //all
        void read(std::istream& is1){
            std::istringstream is_com(read_xyz(is1,l_bohr));
            read_cellp1(is_com,l_bohr);
            afrac();
        }
    };
    
}
using namespace Cell_space;
#endif
