#ifndef CELL_H
#define CELL_H

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <map>

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
    
    virtual std::istream& read(std::istream& is){std::cerr << "read not implement"; return is;}
    virtual std::istream& read_box(std::istream& is){std::cerr << "read_box not implement"; return is;}
    virtual std::istream& read_atoms(std::istream& is){std::cerr << "read_atoms not implement"; return is;}
    virtual std::istream& read_wans(std::istream& is){std::cerr << "read_wans not implement"; return is;}
    virtual std::ostream& write(std::ostream& os){std::cerr << "write not implement"; return os;}
    virtual std::ostream& write_box(std::ostream& os){std::cerr << "write not implement"; return os;}
    virtual std::ostream& write_atoms(std::ostream& os){std::cerr << "write not implement"; return os;}
    virtual std::ostream& write_wans(std::ostream& os){std::cerr << "write not implement"; return os;}
};

#endif
