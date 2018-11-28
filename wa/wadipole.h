#ifndef WADIPOLE_H
#define WADIPOLE_H

#include "andipole.h"
#include "wahbond.h"

class Wadipoles: public Wahbonds{

public:
    typedef typename Cell::position_type position_type;
    
protected:
    position_type tot_dipole;
    vector<shared_ptr<Dipole> > dipoles; // in e*A/ps
    
public:
    static double OW_distance;
    
    Wadipoles() = default;
    Wadipoles(Cell& in_cel, bool ifsave=0);
    //get information
    int nwannier(int i) const;
    int get_wannier(int i,int j) const;
    position_type get_tot_dipole() const;
    position_type get_dipole(int i) const;
    
    void write_each_wannier(ostream& os) const;
    void write_each_dipole(ostream& os) const;
    
    //set
    void set_dipole(int i,position_type D, int in_type=0);
    
    void init_dipoles();
    void set_dipoles_wannier();
    /* type
        0 total dipoles
     
     */
    void set_dipoles_dipole(int type=0);
    void unset_dipoles();
    
private:
    //test
    bool allocate_dipoles(size_t i=0) const;
    void test(bool condition,string information) const;
    //set
    void set_dipoles_total_dipole();
};

class WadipolesT{

public:
    typedef typename Wadipoles::position_type position_type;
    
protected:
    double time = -1;
    double dt   = -1;
    //vector<double> time;
    vector<position_type> dipole,vdipole;
    
public:

    static void Routine();
    
    WadipolesT() = default;
    //get
    double get_dt() const;
    shared_ptr<vector<position_type> > get_dipole() const;
    shared_ptr<vector<position_type> > get_vdipole() const;
    void write_vdipole(ostream& os) const;

    //set
    void clear();
    void add(const Wadipoles& WD);
    void add(const position_type& D, double in_time);
    void read_vdipole(istream& is);
    
    //operations
    void calculate_vdipole();

protected:
    //set
    void set_time(double nt);
    
    //test
    //bool allocate_time(size_t i=0) const;
    bool allocate_dt() const;
    bool allocate_dipole(size_t i=0) const;
    bool allocate_vdipole(size_t i=0) const;
    void test(bool condition,string information) const;
    void log(string information) const;
};

inline int Wadipoles::nwannier(int i) const{
#ifndef N_TEST
    test(allocate_dipoles(i),"ndipole::out_of_range");
#endif
    return dipoles[i]->n();
}
inline int Wadipoles::get_wannier(int i,int j) const{
#ifndef N_TEST
    test(allocate_dipoles(i),"get_index_wannier::out_of_range");
#endif
    return dipoles[i]->get_wannier(j);
}
inline Wadipoles::position_type Wadipoles::get_tot_dipole() const{
    return tot_dipole;
}
inline Wadipoles::position_type Wadipoles::get_dipole(int i) const{
#ifndef N_TEST
    test(allocate_dipoles(i),"get_dipole::out_of_range");
#endif
    return dipoles[i]->get_dipole();
}
inline void Wadipoles::set_dipole(int i,position_type D, int in_type){
#ifndef N_TEST
    test(allocate_dipoles(i),"set_dipole::out_of_range");
#endif
    dipoles[i]->set_dipole(D,in_type);
}
inline void Wadipoles::unset_dipoles(){
    dipoles.clear();
}
inline bool Wadipoles::allocate_dipoles(size_t i) const{
    return i<dipoles.size() and dipoles[i];
}
//WadipolesT::
inline double WadipolesT::get_dt() const{
#ifndef N_TEST
    test(allocate_dt(),"get_dt::no_dt");
#endif
    return dt;
}
inline shared_ptr<vector<typename WadipolesT::position_type> > WadipolesT::get_dipole() const{
#ifndef N_TEST
    test(allocate_dipole(),"get_dipole::no_dipole");
#endif
    return make_shared<vector<position_type> >(dipole);
}
inline shared_ptr<vector<typename WadipolesT::position_type> > WadipolesT::get_vdipole() const{
#ifndef N_TEST
    test(allocate_vdipole(),"get_vdipole::no_vdipole");
#endif
    return make_shared<vector<position_type> >(vdipole);
}
inline void WadipolesT::clear(){
    //time.clear();
    dipole.clear();
    vdipole.clear();
}
//inline bool WadipolesT::allocate_time(size_t i) const{
//    return i<time.size();
//}
inline bool WadipolesT::allocate_dt() const{
    return dt!=-1;
}
inline bool WadipolesT::allocate_dipole(size_t i) const{
    return i<dipole.size();
}
inline bool WadipolesT::allocate_vdipole(size_t i) const{
    return i<vdipole.size();
}

#endif
