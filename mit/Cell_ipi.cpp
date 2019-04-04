#include "Cell_ipi.h"

void Cell_ipi_c::read(std::istream& is){
    switch (cell_type) {
        case 1:
            read_type1(is);
            break;
        default:
            std::cerr << "no read implemented for Cell_ipi_c::type " << cell_type << "! Quit!"; std::exit(0);
            break;
    }
}
void Cell_ipi_c::skip(std::istream& is) const {
    switch (cell_type) {
        case 1:
            skip_type1(is);
            break;
        default:
            std::cerr << "no skip implemented for Cell_ipi_c::type " << cell_type << "! Quit!"; std::exit(0);
            break;
    }
}

void Cell_ipi_c::read_type1(std::istream& is){
    std::string tmps;
    double a,b,c;
    is >> na;
    is >> tmps >> tmps >> a >> b >> c >> tmps >> tmps >> tmps >> tmps >> tmps >> tmps >> c_snapshot;
    is.ignore(1000,'\n');
    read_box(a,b,c);
    init_positions(na);
    for(int i=0; i!=na; ++i){
        is >> tmps >> a >> b >> c;
        types.push_back(tmps);
        read_positions(i,a,b,c);
    }
    set_unit(l_bohr);
}
void Cell_ipi_c::skip_type1(std::istream& is) const{
    double line_count,ipi_unit=1;
    is >> line_count;
    for(int i=0; i!=line_count+2; ++i)
        is.ignore(1000,'\n');
}
