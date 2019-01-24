#ifndef Molecule_h
#define Molecule_h

#include <memory>
#include <vector>

#include "Cell_position.h"

class molecule{
    
    std::vector<std::shared_ptr<position> > atoms_ptrv;
    std::vector<std::shared_ptr<position> > wans_ptrv;
    
public:
    molecule() = default;
    molecule(const molecule&) = default;
    molecule(molecule&&) = default;
    
    std::vector<std::shared_ptr<position> >& atoms(){return atoms_ptrv;}
    std::vector<std::shared_ptr<position> >& wans(){return wans_ptrv;}
};

#endif

