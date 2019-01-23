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
    
public:
    std::shared_ptr<box> box_ptr;
    std::vector<std::shared_ptr<position> > atoms_ptrv;
    std::vector<std::shared_ptr<position> > wans_ptrv;
    std::map<std::string,std::vector<std::shared_ptr<molecule> > > mols_ptrv;
    
    virtual void set_box(double a,double b,double c){
        box_ptr = std::make_shared<box>(Matrix3<double>(a,0,0,0,b,0,0,0,c));
    }
    
    virtual istream& read(istream&){std::cerr << "read not implement"};
    virtual istream& read_box(istream&){std::cerr << "read_box not implement"};
    virtual istream& read_atoms(istream&){std::cerr << "read_atoms not implement"};
    virtual istream& read_wans(istream&){std::cerr << "read_wans not implement"};
    virtual ostream& write(istream&){std::cerr << "write not implement"};
    virtual ostream& write_box(istream&){std::cerr << "write not implement"};
    virtual ostream& write_atoms(istream&){std::cerr << "write not implement"};
    virtual ostream& write_wans(istream&){std::cerr << "write not implement"};
}

#endif
