#ifndef MOLECULE_H
#define MOLECULE_H

#include "Cell.h"

class Molecule: public Atom_group{
    
    Molecule(const Cell& in_cel,std::vector<int> in_indexs=std::vector<int>()):Atom_group(in_cel,in_indexs){}
    Molecule(const Molecule&) = default;
    Molecule(Molecule&&) = default;
};

class Molecule_manip: public Cell_manip{
public:
    Molecule_manip(const Cell& in_cel):Cell_manip(in_cel){}
    

};
#endif
