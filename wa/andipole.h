#ifndef ANDIPOLE_H
#define ANDIPOLE_H

#include "ccell.h"
class Dipole;

class Dipole{
  
public:
    typedef typename Cell::position_type position_type;
    
protected:
    vector<int> wannier;
    int dipole_type = -1; // -1 not allocate;
    position_type dipole; // [eA]
    
public:
    Dipole() = default;
    
    //get
    int n() const;
    int get_wannier(int i) const;
    position_type get_dipole() const;
    bool check_wannier(int i) const;
    
    void write_wannier(ostream& os) const;
    void write_dipole(ostream& os,double unitconv = 1) const;
    
    //set
    void clear();
    void push_back_wannier(int i);
    void swap_wannier(int i,int j);
        // set_wannier for cel->apos(type,index) with distance cutoff
    void set_wannier(shared_ptr<Cell> cel, int type,int index, double distance);
    void set_wannier(shared_ptr<Cell> cel, string type,int index, double distance);
        // set_dipole and dipole type
    void set_dipole(position_type D, int in_type=0); 
    
private:
    //test
    void test(bool condition,string information) const;
    bool allocate_wannier(size_t i=0) const;
    bool allocate_dipole() const;
    
};

inline int Dipole::n() const{
    return wannier.size();
}
inline int Dipole::get_wannier(int i) const{
#ifndef N_TEST
    test(allocate_wannier(i),"get_wannier::out_of_range");
#endif
    return wannier.at(i);
}
inline typename Dipole::position_type Dipole::get_dipole() const{
#ifndef N_TEST
    test(allocate_dipole(),"get_dipole::no_dipole_allocate");
#endif
    return dipole;
}
inline void Dipole::clear(){
    wannier.clear();
    dipole_type = -1;
}
inline void Dipole::push_back_wannier(int i){
    wannier.push_back(i);
}
inline void Dipole::set_dipole(position_type D, int in_type){
    dipole = D;
    dipole_type = in_type;
}
inline bool Dipole::allocate_wannier(size_t i) const{
    return i<wannier.size();
}
inline bool Dipole::allocate_dipole() const{
    return dipole_type!=-1;
}
#endif
