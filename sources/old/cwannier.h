#ifndef CWANNIER_H
#define CWANNIER_H

#include "gheader.h"

/*
 WANNIER DATA
 type:
 intialize:
 Wannier(int in_na = -1);
 input:
 output:

*/
class Wannier;
Inputfile& operator>>(Inputfile &,Wannier &);

class Wannier: public Position_data {
    friend Inputfile& operator>>(Inputfile &,Wannier &);
    //friend class Cell;
    
protected:
    //string      name    =   "wannier";      // name of the atom
    //int         na      =   -1;           // numbers of the atom
    //bool        fractional = false;       // whether fractional
    //shared_ptr<vector<string> >  labels =   make_shared<vector<string> >(); // label informations
    //shared_ptr<vector<position_type>    position   = make_shared<vector<Vector3<double> > >();   // position of atoms in unit of A  or fractional
    
public:
    //Wannier() = default;
    Wannier(int in_na = -1);
    
    //get information
    shared_ptr<Wannier> save() const;
    
    //set information
    void skip(Inputfile &inf) const;
    void read_wfc(istream& is);  // assume  na was set and is start at the first wfc line
    void skip_wfc(istream& is) const;
protected:
    
};


#endif
