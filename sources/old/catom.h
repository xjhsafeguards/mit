#ifndef CATOM_H
#define CATOM_H

#include "gheader.h"

/*
 ATOM TYPE:
 
 initialize:
    Atom();
    Atom(string in_name, int in_na=-1, double in_mass=-1, double in_charge=-1);
 
 input:
    Atom& read_pos();
 
 output:
    double get_mass();
    doubel get_charge();
 
 */

class Atom;
Inputfile& operator>>(Inputfile &,Atom &);
Outputfile& operator<<(Outputfile &, const Atom &);

class Atom : public Position_data {
    friend Inputfile& operator>>(Inputfile &,Atom &);
    friend Outputfile& operator<<(Outputfile &, const Atom &);
    //friend class Cell;
    
private:
    //string      name    =   "default";      // name of the atom
    //int         na      =   -1;           // numbers of the atom
    //bool        fractional = false;       // whether fractional
    //shared_ptr<vector<string> >  labels =   make_shared<vector<string> >(); // label informations
    //shared_ptr<vector<position_type>    position   = make_shared<vector<Vector3<double> > >();   // position of atoms in unit of A  or fractional
    
    double      mass    =   -1;           // in unit of C/12
    double      charge  =   -1;           // in unit of e
    shared_ptr<vector<int> >                 index      = make_shared<vector<int> >(); // the index of atom in a big cell
    
public:
    Atom() = default;
    Atom(string in_name, int in_na=-1, double in_mass=-1, double in_charge=-1);
    //get information
    inline double get_mass() const;
    inline double get_charge() const;
    shared_ptr<Atom> save() const;
    
    void write_in(ostream& os) const;
    void write_in(ostream& os, double unitconv) const;
    void write_POSCAR(ostream& os) const;
    void write_POSCAR(ostream& os, double unitconv) const;
    
    //set information
    inline void set_mass(double in_mass);
    inline void set_charge(double in_charge);
    void skip(Inputfile& inf) const;
    void set_label_file(Inputfile& inf);
    
    void read_pos(istream& is);  // assume na was set and is start at the first atom line
    void skip_pos(istream& is) const;
    void read_cif(istream& is);  // assume name and na was set
    void read_xyz(istream& is);  // assume name and na was set
    void skip_xyz(istream& is) const;
    
    
private:
    //test
    inline bool allocate_index(size_t i=0) const;
};

inline double Atom::get_mass() const {
#ifndef N_TEST
    test_init(mass==-1,name+":mass");
#endif
    return mass;
}
inline double Atom::get_charge() const {
#ifndef N_TEST
    test_init(charge==-1,name+":charge");
#endif
    return charge;
}
inline void Atom::set_mass(double in_mass){
    mass = in_mass;
}
inline void Atom::set_charge(double in_charge){
    charge = in_charge;
}
inline bool Atom::allocate_index(size_t i) const {
    return i<index->size();
}

#endif
