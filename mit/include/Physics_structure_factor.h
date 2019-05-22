#ifndef PHYSICS_STRUCTURE_FACTOR_H
#define PHYSICS_STRUCTURE_FACTOR_H

#include <vector>
#include <iterator>
#include <cassert>

#include "Math_distributionfunction.h"

class Structure_factor{
    
    typedef std::vector<double>::const_iterator g_citer;
    
    int N=-1;
    double dq=-1;
    std::vector<double> S_q,_q;
    
    void compute_sturcture_factor(g_citer begin, g_citer end, double dr, std::back_insert_iterator<std::vector<double>> result, double rou=1, double r0=0, int Ndr=-1);
    void compute_q(std::back_insert_iterator<std::vector<double>> result);
    
public:
    
    Structure_factor(){}
    Structure_factor(const Distributionfunction& in_D, double rou=1, double r0=0, int Ndr=-1);
    
    void Solve(g_citer begin, g_citer end, double dr, double rou=1, double r0=0, int Ndr=-1);
    
    const std::vector<double>& get_s() {return S_q;}
    const std::vector<double>& get_q() {return _q;}
};

#endif
