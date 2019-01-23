#include "Cell_Wannier90.h"

istream& cell_wannier90::read(istream&){
    double line_count;
    is >> line_count;
    is.ignore(500,'\n');
    is.ignore(500,'\n');
    std::shared_ptr<position> tmp_atom = std::make_shared<position>
    
}
