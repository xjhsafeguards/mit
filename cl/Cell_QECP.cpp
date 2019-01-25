#include "Cell_QECP.h"

std::istream& cell_qecp::read_box(std::istream& is){
    is >> snapshot >> time;
    box_ptr = std::make_shared<box>(is);
    box_ptr->set_unit(l_bohr);
    return is;
}

std::istream& cell_qecp::read_atoms(std::istream& is){
    assert(box_ptr);
    double tmps,tmpt;
    is >> tmps >> tmpt;
    assert(snapshot==tmps);
    assert(time==tmpt);
    for(int i=0;i<na.size();++i){
        for(int j=0; j!=na[i]; ++j){
            std::shared_ptr<position> tmp = std::make_shared<cart_position>();
            tmp->set_box_ptr(box_ptr);
            tmp->set_type(types[i]);
            tmp->set_position(is);
            tmp->set_unit(l_bohr);
            atoms_ptrv.push_back(tmp);
        }
    }
    //std::cout << "check0: " << atoms_ptrv.size() << std::endl;
    return is;
}
