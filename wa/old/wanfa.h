//Analyzing wannierfile
#ifndef WANFA_H
#define WANFA_H

#include "wannierfile.h"
#include "input.h"

class Wanfa;
class Wanfat;

class Wanfa
{
    friend class Wanfat;
    
public:
    Wanfa(const Wannierfile& WF_in, const Cellfile& Cel_in);
    Wanfa(const Wannierfile& WF_in, const Cellfile& Cel_in, const Input& INPUT_in);
    
    void Read_dipole();
    ostream& Print_dipole(ostream& os) const;
    
private:
    const Wannierfile& WF;
    const Cellfile& Cel;

    double T;     // Temperature K
    double V;     // Volume m^3
    
    Vector3<double> dipole; //debye
    

    //double V = Cel.Show_cell_volume().pow(INPUT.unitconv,3)* 1e-30; // Volume m^3
    
    //double uc = 36.0*PI/1.38/9.0 *0.20819434*0.20819434*1.6*1.6/2.418884326505/2.418884326505*pow(10.0,13)/INPUT.T/INPUT.Celldm.x/INPUT.Celldm.y/INPUT.Celldm.z/pow(INPUT.unitconv,3);
    //double uc = Miu*C/(3*Kb*SIT*SIV)*D2C*D2C/A2S/A2S*P2S*M2C;
};

class Wanfat
{
public:
    
    Wfat();
    Wfat(const Input& INPUT_in);
private:
    
}

#endif
