#ifndef MOLECULE_H
#define MOLECULE_H

#include <memory> // shared_ptr make_shared
#include <vector> // vector
#include "Cell.h"

class Molecule;
class Molecule_group;

class Molecule: public Atom_group{
    friend class Molecule_group;
public:
    Molecule(const Cell& in_cel,std::vector<int> in_indexs=std::vector<int>()):Atom_group(in_cel,in_indexs){}
    Molecule(const Atom_group& in_AG):Atom_group(in_AG){}
    
    Molecule(const Molecule&) = default;
    Molecule(Molecule&&) = default;
    //using Atom_group::Atom_group;
    
};

class Molecule_group{
public:
    typedef std::shared_ptr<Molecule> molecule_type;
    typedef std::vector<std::shared_ptr<Molecule> > data_type;
protected:
    const Cell& cel;
    data_type molecules;
public:

    Molecule_group(const Cell& in_cel, data_type in_mols=data_type()):cel(in_cel),molecules(in_mols){}
    Molecule_group(const Molecule_group&) = default;
    Molecule_group(Molecule_group&&) = default;
    
    //add new molecules
    void add(const Molecule& m) {molecules.push_back(new_molecule(m));}
    void add(molecule_type mp) {molecules.push_back(mp);}
    
    //iterations
    data_type::iterator begin() {return molecules.begin();}
    data_type::iterator end() {return molecules.end();}
    
    //generate corresponding molecules
    virtual void read() {std::cerr << "no read implemented in Molecule_group! Quit!"; std::exit(0);}
protected:
    
    molecule_type new_molecule() {return std::move(std::make_shared<Molecule>(cel));}
    molecule_type new_molecule(const Molecule& m) {return std::move(std::make_shared<Molecule>(m));}
};
#endif
