#include "Molecule_water.h"

double water_parameter::OH_distance = 1.26;
double water_parameter::OH_distance_tol = 2;
double water_parameter::OW_distance = 0.8;

void Water_group::read(){
    for(const auto& pO : cel.atoms("O")){
        molecule_type tmp_ptr = new_molecule();
        tmp_ptr->add(*pO);
        for(const auto& pH: cel.atoms("H"))
            if(pO->distance(*pH) < water_parameter::OH_distance)
                tmp_ptr->add(*pH);
        //if(tmp_ptr->size()>3)
        std::sort(tmp_ptr->begin()+1,tmp_ptr->end(),[&, pO](std::shared_ptr<Atom>&h1,std::shared_ptr<Atom>&h2){return pO->distance(*h1)<pO->distance(*h2);});
        add(tmp_ptr);
    }
}
