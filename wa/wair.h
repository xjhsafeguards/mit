#ifndef WAIR_H
#define WAIR_H

#include "wadipole.h"
#include "mIR.h"

class WaIR: public Basic_Class{
public:
    static void Routine();
    
private:
    static void C0();  // calculate general IR
    static void C80(); // calculate intermolecule IR with water distance between parameter[0,1]
    static void C81(); // calculate C80 plus INPUT.index Hbond num to first shell
    static void C82(); // calculate C80 plus INPUT.index Hbond num to center
    static void C90(); // calculate intermolecule IR with water has index hb
    static void TEST();

    //test
    static bool in_range(double d);    // test if d between parameter[0,1];
    //standard
};

inline bool WaIR::in_range(double d){
    if(INPUT.parameter.size()<2)
        throw(runtime_error("WaIR::in_range::no_input_parameter"));
    //cout << INPUT.parameter[0] << " " << INPUT.parameter[1] << endl;
    return (INPUT.parameter[0]<d and d<INPUT.parameter[1]);
}

#endif
