#ifndef GPHY_H
#define GPHY_H

#include "gfun.h"

class Unit
//Value with units
{
public:
    static constexpr double Kb = 1.38064852*1e-23;  // boltzman constant J/K
    static constexpr double C  = 299792458;      // speed of light m/s
    static constexpr double Miu = 4*3.141592653589793238460*1e-7;    // magnetic constant H/m or kg/s^2/A^2
    static constexpr double Debye = 0.20819434*1.6*1e-19*1e-10;  // unit conv of dipole from debye to C*M
    static constexpr double AU_t = 2.418884326505*1e-17; // unit conv of time from a.u. time to s
    static constexpr double Ps = 1e-12;     // unit conv of time from ps to s
    static constexpr double M2C = 1e-5;      // unit conv of final output from m^-1 to 10^3 cm^-1
    static constexpr double Bohr =  5.2917721067e-11; // unit conv of Bohr to m
    static constexpr double Bohr2A = 0.52917721067;     // unit conv of Bohr to Angstrom
    static constexpr double Cm_12Ps_1 = 6*3.141592653589793238460/100; // unit conv of cm-1 to 2PIps-1
    static constexpr double EA = 1.6*1e-19*1e-10;    // unit conv of eA to C*M
    static constexpr double A3 = 1e-30;
    static double unitconv(string unit);
};


#endif
