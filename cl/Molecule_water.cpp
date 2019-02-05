#include "Molecule_water.h"

double water_parameter::OH_distance = 1.26;
double water_parameter::OW_distance = 0.8;

void water_manip::read_atoms(const std::vector<std::shared_ptr<position> >& atoms_ptrv,std::vector<std::shared_ptr<molecule> >& out_moleculev){
    assert(atoms_ptrv.size());
    for(const auto& aptr1 : atoms_ptrv){
        if(aptr1->check_type("O")){
            std::shared_ptr<molecule> tmp_ptr = std::make_shared<molecule>();
            tmp_ptr->atoms().push_back(aptr1);
            for(const auto& aptr2 : atoms_ptrv){
                if(aptr2->check_type("H") and aptr1->distance(*aptr2) < water_parameter::OH_distance){
                    tmp_ptr->atoms().push_back(aptr2);
                }
            }
            // if H2O has more than 2 H store H in distance order
            if(tmp_ptr->atoms().size()>3){
                auto atom_o = tmp_ptr->atoms()[0];
                sort(tmp_ptr->atoms().begin()+1,tmp_ptr->atoms().end(),[&, atom_o](std::shared_ptr<position>&h1,std::shared_ptr<position>&h2){return atom_o->distance(*h1)<atom_o->distance(*h2);});
                // H-O-H angle should be bigger than 80
                //assert(atom_o->angle(*(tmp_ptr->atoms()[1]),*(tmp_ptr->atoms()[2]))>80);
                //std::cout << tmp_ptr->atoms().size()-1 << "H molecules found: " << atom_o->angle(*(tmp_ptr->atoms()[1]),*(tmp_ptr->atoms()[2])) << '\n';
            }
            out_moleculev.push_back(tmp_ptr);
        }
    }
}

void water_manip::read_atoms(cell& in_cell){
    read_atoms(in_cell.atoms(),in_cell.mols("H2O"));
}

void water_manip::read_wans(const wans_type& wans_ptrv,mols_type& out_moleculev){
    assert(wans_ptrv.size());
    assert(out_moleculev.size());
    for(const auto& waterp : out_moleculev){
        for(const auto& wanp : wans_ptrv){
            if( waterp->atoms()[0]->distance(*wanp) < water_parameter::OW_distance){
                waterp->wans().push_back(wanp);
            }
        }
        if( waterp->wans().size() != 4 )
        {
            std::cerr << "Something went wrong while reading wannier center in water molecule: " << waterp->atoms()[0]->cart() <<std::endl;
            std::cerr << "It has " << waterp->wans().size() << " wannier centers (should be 4)" << std::endl;
            std::cerr << "Try to change a wannier center lenth cutoff!"<< std::endl;
            exit(1);
        }
    }
}
void water_manip::read_wans(cell& in_cell){
    read_wans(in_cell.wans(),in_cell.mols("H2O"));
}
void water_manip::read(cell& in_cell){
    read_atoms(in_cell);
    //std::cout << "check 4" << std::endl;
    if(in_cell.wans().size())
        read_wans(in_cell);
}

void water_manip::sort_wans(mols_type& out_moleculev){
    assert(out_moleculev.size());
    assert(out_moleculev[0]->wans().size());
    for(const auto& waterp : out_moleculev){
        const auto& atom_o = waterp->atoms()[0];
        //sort wans in the descending order of distance to O
        sort(waterp->wans().begin(),waterp->wans().end(),[&, atom_o](std::shared_ptr<position>&w1,std::shared_ptr<position>&w2){return atom_o->distance(*w1)>atom_o->distance(*w2);});
        //if water has more than 1 H, sort wans in the order of H
        if(waterp->atoms().size()>2 and waterp->atoms()[1]->distance(*(waterp->wans()[0]))>waterp->atoms()[1]->distance(*(waterp->wans()[1])))
            swap(waterp->wans()[0],waterp->wans()[1]);
    }
}

void water_manip::sort_wans(cell& in_cell){
    sort_wans(in_cell.mols("H2O"));
}
