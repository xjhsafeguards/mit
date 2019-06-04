#include "Cell.h"
#include "Physics_unit.h"

void Cell_qecp::read_cellp(std::istream& is){
    //assert(is.good());
    is >> snapshot >> time;
    cellp.load(is);
    cellp.set_unit(l_bohr);
}
void Cell_qecp::skip_cellp(std::istream& is) const{
    //assert(is.good());
    for(int i=0;i!=4;++i)
        is.ignore(1000,'\n');
}
void Cell_qecp::read_atoms(std::istream& is){
    //assert(is.good());
    is >> snapshot >> time;
    assert(na!=-1);
    atoms.init(na);
    for(int i=0; i!=na; ++i){
        atoms.load1(i,is);
    }
    atoms.set_unit(l_bohr);
    atoms.frac(cellp);
}
void Cell_qecp::skip_atoms(std::istream& is) const{
    //assert(is.good());
    for(int i=0;i!=na+1;++i)
        is.ignore(1000,'\n');
}
void Cell_qecp::read(std::istream& is1,std::istream& is2){
    read_cellp(is1);
    read_atoms(is2);
}
void Cell_qecp::skip(std::istream& is1,std::istream& is2) const {
    skip_cellp(is1);
    skip_atoms(is2);
}
