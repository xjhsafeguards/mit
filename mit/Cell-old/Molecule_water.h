#ifndef Molecule_water_h
#define Molecule_water_h

#include <algorithm>
#include <vector>

#include "Molecule.h"
#include "Cell.h"
#include "Moleculemanip.h"

namespace water_parameter {
    extern double OH_distance;
    extern double OH_distance_tol;
    extern double OW_distance;
    }

// Water structure:
//      O H H ...
//      H1 H2 lone lone

class water_manip: public molecule_manip{
public:
    virtual void read_atoms(const atoms_type&,mols_type& out_moleculev);
    virtual void read_atoms(cell& in_cell);
    virtual void read_wans(const wans_type&,mols_type& out_moleculev);
    virtual void read_wans(cell& in_cell);
    virtual void read(cell& in_cell);
    
    virtual void sort_wans(mols_type& out_moleculev);
    virtual void sort_wans(cell& in_cell);
};

#endif
