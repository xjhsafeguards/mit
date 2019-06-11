#include "Cell.h"


void Cell::read_cel(std::istream& is){
    //assert(is.good());
    is >> snapshot >> time;
    cellp.load(is);
    cellp.set_unit(l_bohr);
}
void Cell::skip_cel(std::istream& is) const{
    //assert(is.good());
    for(int i=0;i!=4;++i)
        is.ignore(1000,'\n');
}
void Cell::read_pos(std::istream& is){
    //assert(is.good());
    is >> snapshot >> time;
    assert(na!=-1);
    atoms.init(na);
    for(int i=0; i!=na; ++i){
        atoms.load1(i,is);
    }
    atoms.set_unit(l_bohr);
}
void Cell::skip_pos(std::istream& is) const{
    //assert(is.good());
    for(int i=0;i!=na+1;++i)
        is.ignore(1000,'\n');
}
std::string Cell::read_xyz(std::istream& is,double in_unit){
    std::string tmps,result;
    is >> na;
    is.ignore(1000,'\n');
    getline(is,result);
    atoms.init(na);
    for(int i=0; i!=na; ++i){
        is >> tmps;
        types.push_back(tmps);
        atoms.load1(i,is);
    }
    atoms.set_unit(in_unit);
    return result;
}
void Cell::skip_xyz(std::istream& is) const{
    double line_count;
    is >> line_count;
    for(int i=0; i!=line_count+2; ++i)
        is.ignore(1000,'\n');
}


/*****************************   CELL_IPI    **********************************/
void Cell_ipi::read_cellp1(std::istream& is,double in_unit){
    std::string tmps;
    double a,b,c;
    is >> tmps >> tmps >> a >> b >> c >> tmps >> tmps >> tmps >> tmps >> tmps >> tmps >> snapshot;
    cellp.load(a,b,c);
    cellp.set_unit(in_unit);
}
