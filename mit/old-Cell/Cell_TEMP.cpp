#include "Cell_TEMP.h"

std::istream& cell_temp1::read_atoms(std::istream& is){
    assert(box_ptr);
    std::string tmpstring;
    is >> tmpstring >> snapshot;
    for(int i=0; i!=190; ++i){
        std::shared_ptr<position> tmp = std::make_shared<cart_position>();
        tmp->set_box_ptr(box_ptr);
        tmp->set_type(is);
        tmp->set_position(is);
        is.ignore(500,'\n');
        atoms_ptrv.push_back(tmp);
    }
    //std::cout << "check0: " << atoms_ptrv.size() << std::endl;
    return is;
}
