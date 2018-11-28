#ifndef DQEPOS_H
#define DQEPOS_H

#include "gheader.h"

class DQEpos{
    typedef vector<double>  Positions;
    typedef double          Time;
public:
    // I/O
    
private:
    Positions   positions;
    Time        time;
};

#endif
