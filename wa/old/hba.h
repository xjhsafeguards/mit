// H bond analyze
#ifndef HBA_H
#define HBA_H

#include "H_bondfile.h"

class Hba
{
public:
    
    Hba(Hbondfile& HBF_in, Cellfile& Cel_in);
    ~Hba();
    
    static void Routine(); //1 rotational correlation time (function)
    
    /* basic values */
    Hbondfile& HBF; //the H_bond file to analyze
    Cellfile& Cel; //the corresponding cellfile
    int n_HB; //total number of HBs in the system(only count for donate)
    Vector3<double> rOH; //the average vector from O to H in Hbond (SI A)
    
    /* basic functions */
    void Read_rOH(); // Read rOH from HBF and Cel
    void Print_rOH(ofstream& ofs); // Print rOH into ofs;
    
private:
    
    
};

class Hbafile
//store rOH data over time
{
public:
    
    Hbafile(int N);
    ~Hbafile();
    
    double dt; //time interval of snapshot
    int N; //numbers of snapshot contains start from 1
    
    Vector3<double> *rOH; //rOH for different time
    
    /* rOH time correlation function */
    bool allocate_rOH_tcf;
    int n_rOH_tcf; // number of tcf calculated do not count <0*0>
    int n_jump; // dt_tcf = dt * n_jump
    double *rOH_tcf; //<P2(rOH(t)*rOH(0))>/<P2(rOH(0)^2)> and the first one is <P2(rOH(0)^2)>
    void Calculate_rOH_tcf(int n); // calculate rOH_tcf with max n step
    void Calculate_rOH_tcf(int n, ofstream& ofs);
    double t2; //rotational correlation time (SI ps)
    void Calculate_t2(ofstream& ofs);
    
    /* basic functions */
    static double Time_correlation_Legendre2(Vector3<double> *y, int n_max, int n_step);
    
};

#endif
