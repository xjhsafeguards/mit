#ifndef PCONST_H
#define PCONST_H

#include "punit.h"

namespace physcis_const {
    const double Kb = 1.38064852*1e-23*physics_unit::e_j/physics_unit::t_k;  // boltzman constant ev/K (J/K)
    const double C  = 299792458*physics_unit::l_m/physics_unit::t_s;        // speed of light angstrom/ps (m/s)
    const double Miu = 4*3.141592653589793238460*1e-7;    // magnetic constant H/m or kg/s^2/A^2
    const double Na = 6.02214082e23; // avogadro constant
}

using physcis_const::Kb;
using physcis_const::C;
using physcis_const::Miu;
using physcis_const::Na;


#endif
