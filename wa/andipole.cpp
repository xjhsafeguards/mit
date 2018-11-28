#include "andipole.h"

bool Dipole::check_wannier(int i) const{
    auto result = find(wannier.cbegin(),wannier.cend(),i);
    return result!=wannier.cend();
}
void Dipole::write_wannier(ostream& os) const{
    for(const auto& wi: wannier)
        os << setw(5) << wi;
}
void Dipole::write_dipole(ostream& os,double unitconv) const{
    os << dipole/unitconv;
}
void Dipole::swap_wannier(int i,int j){
#ifndef N_TEST
    test(allocate_wannier(i) and allocate_wannier(j),"swap_wannier:out_of_range");
#endif
    int tmp = wannier[i];
    wannier[i] = wannier[j];
    wannier[j] = tmp;
}
void Dipole::set_wannier(shared_ptr<Cell> cel, int type, int index, double distance){
    wannier.clear();
    for(int j=0; j<cel->wnum(); j++)
        if( cel->distance(type,index,j) < distance)
            wannier.push_back(j);
}
void Dipole::set_wannier(shared_ptr<Cell> cel, string type, int index, double distance){
    wannier.clear();
    for(int j=0; j<cel->wnum(); j++)
        if( cel->distance(type,index,j) < distance)
            wannier.push_back(j);
}
void Dipole::test(bool condition,string information) const{
    if(!condition)
        throw(runtime_error("Dipole::Fail_"+information));
}
