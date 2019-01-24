#include "Cell_QECP.h"

std::istream& cell_qecp::read_atoms(std::istream& is){
    assert(box_ptr);
    is >> snapshot >> time;
/*
    for(int i=0; i!=line_count; ++i){
        std::shared_ptr<position> tmp = std::make_shared<cart_position>();
        tmp->set_box_ptr(box_ptr);
        tmp->set_type(is);
        tmp->set_position(is);
        if(tmp->check_type("X"))
            wans_ptrv.push_back(tmp);
        else{
            //std::cout << tmp->cart() << std::endl;
            atoms_ptrv.push_back(tmp);
        }
    }
    //std::cout << "check0: " << atoms_ptrv.size() << std::endl;
 */
    return is;
}
