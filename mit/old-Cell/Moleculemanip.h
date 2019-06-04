#ifndef MOLECULE_MANIP_H
#define MOLECULE_MANIP_H

#include "Molecule.h"
#include "Cell.h"

class molecule_manip{
public:
    typedef std::vector<std::shared_ptr<position> > atoms_type;
    typedef std::vector<std::shared_ptr<position> > wans_type;
    typedef std::vector<std::shared_ptr<molecule> > mols_type;
    
    virtual void read_atoms(const atoms_type&,mols_type& out_moleculev){}
    virtual void read_atoms(cell& in_cell){}
    virtual void read_wans(const wans_type&,mols_type& out_moleculev){}
    virtual void read_wans(cell& in_cell){}
    virtual void read(cell& in_cell){}
    
    virtual void sort_wans(mols_type& out_moleculev){}
    virtual void sort_wans(cell& in_cell){}
};

#endif
