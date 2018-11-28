#include "wanfa.h"

Wanfa::Wanfa(const Wannierfile& WF_in, const Cellfile& Cel_in):
    WF(WF_in), Cel(Cel_in)
{}

Wanfa::Wanfa(const Wannierfile& WF_in, const Cellfile& Cel_in,const Input& INPUT_in):
    WF(WF_in), Cel(Cel_in), T(INPUT_in.T),
{}
