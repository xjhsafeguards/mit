#include "Cell_Lammps.h"
#include "Physics_unit.h"

std::istream& cell_lammps::read(std::istream& is){
    assert(is.good());
    std::string tmps;
    switch (cell_type) {
        case 0:
            _read_timestep(is);
            _read_n(is);
            _read_box0(is);
            _read_pos(is);
            break;
        default:
            std::cerr << "Cell_lammps cell_typ " << cell_type << " not implemented" << std::endl;
        break;
    }
    return is;
}

std::istream& cell_lammps::_read_timestep(std::istream& is){
    std::string tmps;
    do{
        std::getline(is,tmps);
    }while(tmps.find("ITEM: TIMESTEP")!=std::string::npos);
    is >> snapshot;
    return is;
}
std::istream& cell_lammps::_read_n(std::istream& is){
    std::string tmps;
    do{
        std::getline(is,tmps);
    }while(tmps.find("ITEM: NUMBER OF ATOMS")!=std::string::npos);
    is >> na;
    return is;
}
std::istream& cell_lammps::_read_box0(std::istream& is){
    std::string tmps;
    do{
        std::getline(is,tmps);
    }while(tmps.find("ITEM: BOX")!=std::string::npos);
    double a,b,c;
    is >> tmps >> a >> tmps >> tmps >> b >> tmps >> tmps >> c;
    is.ignore(1000,'\n');
    set_box(a,b,c);
    box_ptr->set_unit(cell_unit);
    return is;
}
std::istream& cell_lammps::_read_pos(std::istream& is){
    assert(na!=0);
    std::string tmps;
    do{
        std::getline(is,tmps);
    }while(tmps.find("ITEM: ATOMS")!=std::string::npos);
    int itype;
    for(int i=0; i!=na; ++i){
        std::shared_ptr<position> tmp = std::make_shared<cart_position>();
        tmp->set_box_ptr(box_ptr);
        is >> tmps >> itype;
        tmp->set_type(types[itype-1]);
        tmp->set_position(is);
        tmp->set_unit(cell_unit);
        atoms_ptrv.push_back(tmp);
    }
    return is;
}
