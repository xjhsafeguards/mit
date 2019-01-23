#ifndef Molecule_water_h
#define Molecule_water_h


#include <vector>

#include "Molecule.h"

namespace water_parameter {
    double OH_distance = 1.26;
    }

extern void Read_water(const std::vector<std::shared_ptr<position> >&,std::vector<std::shared_ptr<molecule> >&);

#endif
