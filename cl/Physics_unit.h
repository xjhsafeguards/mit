#ifndef PUNIT_H
#define PUNIT_H

namespace physics_unit {
    // time 1 = ps
    const double t_ps = 1;    // pecosecond
    const double t_s = 1e12;   // second
    const double t_au = 2.418884326505*1e5; // atomic unit time
    // length 1 = Angstrom
    const double l_angs = 1;   // angstrom
    const double l_bohr = 0.52917721067; //bhor
    const double l_m = 1e10;    // meter
    const double l_cm = 1e8;    // centimeter
    const double l_mm = 1e7;    // milimeter
    const double l_nm = 10;     // nanometer
    // mass 1 = atomic mass
    const double m_au = 1;
    const double m_g = 6.02214082e23; // gram
    const double m_kg = 1e3*m_g;  // kilogram
    // charge 1 = e
    const double c_e = 1;
    const double c_c = 1/1.60217662e-19;
    // energy 1 = ev
    const double e_ev = 1;      // eV
    const double e_j = 1/1.60217662e-19; // joules
    // dipole 1 = ea
    const double d_ea = 1;      // e*angstrom
    const double d_debye = 0.20819434;  // debye
    const double d_cm = c_c*l_m;      // c*m
    // temperature 1 = K
    const double t_k = 1; // kelven
}

using physics_unit::t_ps;
using physics_unit::t_s;
using physics_unit::t_au;
using physics_unit::l_angs;
using physics_unit::l_bohr;
using physics_unit::l_m;
using physics_unit::l_cm;
using physics_unit::l_mm;
using physics_unit::l_nm;
using physics_unit::m_kg;
using physics_unit::m_g;
using physics_unit::c_e;
using physics_unit::c_c;
using physics_unit::e_ev;
using physics_unit::e_j;
using physics_unit::d_ea;
using physics_unit::d_debye;
using physics_unit::d_cm;
using physics_unit::t_k;


#endif
