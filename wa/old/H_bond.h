#ifndef H_BOND_H
#define H_BOND_H

#include "water.h"

class Hbond
{
    public:
    
    Hbond();
    ~Hbond();
    
    int accept_count;
    int donate_count;
    int HB_count;
    
    vector<int> index_accept;
    vector<int> index_donate;
    vector<double> donate_angle;
    
    public:
    
    void Print(ofstream &ofs);
    //print the a d ia*8 id*4
};

#endif
