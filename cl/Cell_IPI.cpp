#include "Cell_IPI.h"

std::istream& cell_ipi::read(std::istream& is){
    double line_count,ipi_unit=1;
    std::string tmps;
    is >> line_count;
    switch (cell_type) {
        case 0:{
            double a,b,c;
            is >> tmps >> tmps >> a >> b >> c >> tmps >> tmps >> tmps >> tmps >> tmps >> tmps >> snapshot;
            is.ignore(1000,'\n');
            set_box(a,b,c);
            box_ptr->set_unit(l_bohr);
            ipi_unit = l_bohr;
        }
            break;
        default:
            std::cerr << "Cell_ipi cell_typ " << cell_type << " not implemented" << std::endl;
            break;
    }
    for(int i=0; i!=line_count; ++i){
        std::shared_ptr<position> tmp = std::make_shared<cart_position>();
        tmp->set_box_ptr(box_ptr);
        tmp->set_type(is);
        tmp->set_position(is);
        tmp->set_unit(ipi_unit);
        atoms_ptrv.push_back(tmp);
        }
    return is;
}

std::ostream& cell_ipi::write(std::ostream& os){
    
}
