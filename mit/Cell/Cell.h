#ifndef CELL_H
#define CELL_H

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <unordered_map>
#include <utility> // pair
#include <cstdint>

#include "Cell_position.h"
#include "Molecule.h"

class cell;

class cell{
    
protected:
    double time = -1;
    int snapshot = -1;
    std::shared_ptr<box> box_ptr;
    std::vector<std::shared_ptr<position> > atoms_ptrv;
    std::vector<std::shared_ptr<position> > wans_ptrv;
    
    std::map<std::string,std::vector<std::shared_ptr<molecule> > > mols_ptrv;
    std::unordered_map<std::string,std::vector<std::shared_ptr<position> > > _type_atoms_ptrv;
    std::unordered_map<std::uintptr_t,double> _atoms_dis;

public:
    
    std::shared_ptr<box>& boxp(){return box_ptr;}
    std::vector<std::shared_ptr<position> >& atoms(){return atoms_ptrv;}
    std::vector<std::shared_ptr<position> >& wans(){return wans_ptrv;}
    std::map<std::string,std::vector<std::shared_ptr<molecule> > >& mols(){return mols_ptrv;}
    std::vector<std::shared_ptr<molecule> >& mols(std::string mol_name){return mols_ptrv[mol_name];}
    
    virtual void set_box(double a,double b,double c){
        box_ptr = std::make_shared<box>(Matrix3<double>(a,0,0,0,b,0,0,0,c));
    }
    virtual void clear(){
        atoms_ptrv.clear();
        wans_ptrv.clear();
        mols_ptrv.clear();
    }
    double volume() const{
        return box_ptr->volume();
    }
    int ss() const{
        return snapshot;
    }
    
    //type dependent calculations
    void sort_type(){
        for(const auto& atom : atoms_ptrv){
            _type_atoms_ptrv[atom->type()].push_back(atom);
        }
    }
    std::vector<std::shared_ptr<position> >& atom_type(const std::string& in_type){
        return _type_atoms_ptrv.at(in_type);
    }
    //fast pair distance calculations (slow currently)
    void cal_atoms_pair_dis(){
        for(auto it=atoms_ptrv.cbegin();it!=atoms_ptrv.cend();++it)
            for(auto it2=it;it2!=atoms_ptrv.cend();++it2)
                _insert_pair_dis(*it,*it2,_atoms_dis);
    }
    double atoms_dis(const std::shared_ptr<position>& p1,const std::shared_ptr<position>& p2){
        return _atoms_dis.at(reinterpret_cast<std::uintptr_t>(p1.get())*reinterpret_cast<std::uintptr_t>(p2.get()));
    }
    
    //IO
    
    virtual std::istream& read(std::istream& is){std::cerr << "read not implement"; return is;}
    virtual std::istream& read_box(std::istream& is){std::cerr << "read_box not implement"; return is;}
    virtual std::istream& read_atoms(std::istream& is){std::cerr << "read_atoms not implement"; return is;}
    virtual std::istream& read_wans(std::istream& is){std::cerr << "read_wans not implement"; return is;}
    virtual std::istream& skip(std::istream& is){std::cerr << "read not implement"; return is;}
    virtual std::istream& skip_box(std::istream& is){std::cerr << "read_box not implement"; return is;}
    virtual std::istream& skip_atoms(std::istream& is){std::cerr << "read_atoms not implement"; return is;}
    virtual std::istream& skip_wans(std::istream& is){std::cerr << "read_wans not implement"; return is;}
    virtual std::ostream& write(std::ostream& os){std::cerr << "write not implement"; return os;}
    virtual std::ostream& write_box(std::ostream& os){box_ptr->write(os); return os;}
    virtual std::ostream& write_atoms(std::ostream& os){
        for(const auto& ap : atoms_ptrv){
            ap->write(os) << std::endl;
        }
        return os;}
    virtual std::ostream& write_atoms_frac(std::ostream& os){
        for(const auto& ap : atoms_ptrv){
            ap->write_frac(os) << std::endl;
        }
        return os;}
    virtual std::ostream& write_atoms_cart(std::ostream& os){
        for(const auto& ap : atoms_ptrv){
            ap->write_cart(os) << std::endl;
        }
        return os;}
    virtual std::ostream& write_wans(std::ostream& os){std::cerr << "write not implement"; return os;}
    
private:
    void _insert_pair_dis(const std::shared_ptr<position>& p1,const std::shared_ptr<position>& p2, std::unordered_map<std::uintptr_t,double>& in_map){
        in_map.insert(std::make_pair(reinterpret_cast<std::uintptr_t>(p1.get())*reinterpret_cast<std::uintptr_t>(p2.get()),p1->distance(*p2)));
    }
};

#endif
