#include "Cell_qe.h"

void Cell_qecp_c::read_cell(std::istream& is){
    is >> c_snapshot >> c_time;
    double a,b,c,tmp;
    is >> a >> tmp >> tmp >> tmp >> b >> tmp >> tmp >> tmp  >> c;
    read_box(a,b,c);
    set_box_unit(l_bohr);
}
void Cell_qecp_c::skip_cell(std::istream& is) const{
    for(int i=0;i!=4;++i)
        is.ignore(1000,'\n');
}
void Cell_qecp_c::read_pos(std::istream& is){
    is >> c_snapshot >> c_time;
    double a,b,c;
    init_positions(na);
    for(int i=0; i!=na; ++i){
        is >> a >> b >> c;
        read_positions(i,a,b,c);
    }
    set_positions_unit(l_bohr);
}
void Cell_qecp_c::skip_pos(std::istream& is) const{
    for(int i=0;i!=na+1;++i)
        is.ignore(1000,'\n');
}
void Cell_qecp_c::read(std::istream& is1,std::istream& is2){
    read_cell(is1);
    read_pos(is2);
}
void Cell_qecp_c::skip(std::istream& is1,std::istream& is2) const {
    skip_cell(is1);
    skip_pos(is2);
}
