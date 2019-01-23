#include "Molecule_water.h"

void Read_water(const std::vector<std::shared_ptr<position> >& atoms_ptrv,std::vector<std::shared_ptr<molecule> >& out_moleculev){
    for(const auto& aptr1 : atoms_ptrv){
        if(aptr1->check_type("O")){
            std::shared_ptr<molecule> tmp_ptr = std::make_shared<molecule>();
            tmp_ptr->atoms_ptrv.push_back(aptr1);
            for(const auto& aptr2 : atoms_ptrv){
                if(aptr2->check_type("H") and aptr1->distance(*aptr2) < water_parameter::OH_distance){
                    tmp_ptr->atoms_ptrv.push_back(aptr2);
                }
            }
            out_moleculev.push_back(tmp_ptr);
        }
    }

}
