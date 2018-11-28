#ifndef PHY_H
#define PHY_H

#include "gfun.h"
struct Unit;
class Tcf;

class Unit
//Unitcon
{
public:
    static constexpr double Kb = 1.38064852*1e-23;  // boltzman constant J/K
    static constexpr double C  = 299792458;      // speed of light m/s
    static constexpr double Miu = 4*3.141592653589793238460*1e-7;    // magnetic constant H/m or kg/s^2/A^2
    static constexpr double Debye = 0.20819434*1.6*1e-19*1e-10;  // unit conv of dipole from debye to C*M
    static constexpr double AU_t = 2.418884326505*1e-17; // unit conv of time from a.u. time to s
    static constexpr double Ps = 1e-12;     // unit conv of time from ps to s
    static constexpr double M2C = 1e-5;      // unit conv of final output from m^-1 to 10^3 cm^-1
};

class Tcf
//Time correlation fuction
{
public:
    
    Tcf(vector<double>& data_in, vector<double>& data_out,double dt_in, int n_tcf_in = 0, int n_jump_in = 1) :
        data1(data_in),
        tcf(data_out),
        N(data_in.size()),
        dt(dt_in),
        n_jump(n_jump_in),
        n_tcf(n_tcf_in)   {}
    
    ~Tcf(){};
    
private:
    
    const vector<double>& data1; // hold double data to calculate
    
    unsigned N; // size of data
    double dt; // time interval of snapshot
    
    int n_jump; // dt_tcf = dt * n_jump
    int n_tcf;  // number of tcf calculated not including <0*0>, i.e. last index of rOH_tcf
    vector<double>& tcf;
    bool allocate_tcf = false;
    
    void set_up_calculation();
    
public:
    
    void Calculate_tcf();
    void Calculate_tcf_legendre(int n=1);
    
    vector<double>& Show_tcf() const
    {   return tcf;}
    ostream& Print_tcf(ostream &os = cout) const;
    
};

#endif
