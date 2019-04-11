#ifndef MOLECULE_WATER_H
#define MOLECULE_WATER_H

#include <algorithm>
#include <vector>

#include "Molecule.h"

namespace water_parameter {
    extern double OH_distance;
    extern double OH_distance_tol;
    extern double OW_distance;
}

// Water structure:
//      O H H ...
//      H1 H2 lone lone

class Water_group: public Molecule_group{
public:
    using Molecule_group::Molecule_group;
    
    virtual void read();
};
#endif
