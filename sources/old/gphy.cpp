#include "gphy.h"

double Unit::unitconv(string unit){
    if(unit == "bohr")
        return Unit::Bohr2A;
    else if(unit == "cm-1")
        return Unit::Cm_12Ps_1;
    return 1;
}
