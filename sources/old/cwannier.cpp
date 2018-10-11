#include "cwannier.h"

Wannier::Wannier(int in_na): Position_data("wannier",in_na) {
}

Inputfile& operator>>(Inputfile &inf,Wannier &wan){
    if(inf.type() == "wfc") wan.read_wfc(inf.stream());
    else
        wan.test(0,"WANNIER::operator>>::wrong_inputfile_type:" + inf.type());
    return inf;
}
void Wannier::skip(Inputfile &inf) const{
    if(inf.type() == "wfc") skip_wfc(inf.stream());
    else
        test(0,"WANNIER::operator>>::wrong_inputfile_type:" + inf.type());
}

shared_ptr<Wannier> Wannier::save() const{
    shared_ptr<Wannier> tmpwan = make_shared<Wannier>(*this);
    tmpwan->labels = value_copy(labels);
    tmpwan->position = value_copy(position);
    return tmpwan;
}
void Wannier::read_wfc(istream& is){
    try{
        set_label(is);
        fractional = false;
        set_position(is,Unit::Bohr2A);
    }catch(runtime_error err){
        test(0,"WANNIER:" + name + "::read_pos");
    }
}
void Wannier::skip_wfc(istream& is) const{
    try{
        is.ignore(500,'\n');
        skip_position(is);
    }catch(runtime_error err){
        test(0,"WANNIER:" + name + "::skip_pos");
    }
}


