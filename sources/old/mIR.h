#ifndef MIR_H
#define MIR_H

#include "gheader.h"
class IR;
class IR2;
template<typename Tc>
class IR3;

class IR: public Basic_Class {

public:
    typedef Vector3<double> dipole_type;
    
protected:
    double dt   = -1;  // in A
    shared_ptr<vector<dipole_type> > dipole;    // in e*A
    shared_ptr<vector<dipole_type> > vdipole;   // in e*A/ps
    vector<double> tcf;
    vector<double> ir;
    double dw   = -1;       // in 2Pi*ps-1
    double T = 300;         // in K and A^3
    double V = -1;
    double screen_t = 0;
    double screen_alpha = 0;
    
    int tcf_nmax = -1; // tcf step max
    double ft_kmax; // max frequency for ft in 2PI * ps-1

public:
    IR();
    //get
    double get_dt() const;
    const vector<double>& get_tcf() const;
    
    void write_tcf(ostream& os) const;
    void write_ir(ostream& os) const;
    //set
    void unset_tcf();
    void unset_ir();
    
    void set_screen(double st,double sa);
    void set_dt(double in_dt);
    void set_T(double st);
    void set_V(double V);
    void set_dipole(shared_ptr<vector<dipole_type> > in_d);
    void set_vdipole(shared_ptr<vector<dipole_type> > in_v);
    void set_tcf(const vector<double>& in_tcf);
    void set_up(const IR& in_ir);
    //operations
    void calculate_vdipole();
    void calculate_tcf();
    void calculate_ir();
    
protected:
    //operations
    void c_vdipole(shared_ptr<vector<dipole_type> > d, shared_ptr<vector<dipole_type> > v); // dont check d and v
    void c_tcf_init(shared_ptr<vector<dipole_type> > v1, shared_ptr<vector<dipole_type> > v2, int nmax);
    void c_tcf(shared_ptr<vector<dipole_type> > v1, shared_ptr<vector<dipole_type> > v2, int nmax);
    void gaussian_screen_tcf();
    void unit_convert_ir();
    
    //test
    bool allocate_dt() const;
    bool allocate_dipole(size_t i=0) const;
    bool allocate_vdipole(size_t i=0) const;
    bool allocate_tcf(size_t i=0) const;
    bool allocate_ir(size_t i=0) const;
    bool allocate_V() const;
    
    void test(bool condition,string information) const;
    //void log(string information) const;
};


class IR2: public IR{

protected:
    vector<shared_ptr<vector<dipole_type> > > dipole;    // in e*A
    vector<shared_ptr<vector<dipole_type> > > vdipole;   // in e*A/ps
    
public:
    IR2();
    //get
    
    //set
    void set_dipole(const vector<shared_ptr<vector<dipole_type> > > &in_d);
    void push_back_dipole(shared_ptr<vector<dipole_type> >  in_d);
    void push_back_dipole(const vector<dipole_type>&  in_d);
    void set_vdipole(const vector<shared_ptr<vector<dipole_type> > > &in_v);
    void push_back_vdipole(shared_ptr<vector<dipole_type> >  in_v);
    void push_back_vdipole(const vector<dipole_type>&  in_d);
    //operations
    void calculate_vdipole();
    void calculate_tcf(int i=0, int j=0);

protected:
    bool allocate_dipole(size_t i=0, size_t j=0) const;
    bool allocate_vdipole(size_t i=0, size_t j=0) const;
};

template<typename Tc>
class IR3: public IR2 {
    
    vector<shared_ptr<vector<Tc> > >      condition;   // conditions
    
public:
    IR3();
    //set
    void set_condition(const vector<shared_ptr<vector<Tc> > >& in_c);
    void push_back_condition(shared_ptr<vector<Tc> > in_c);
    void push_back_condition(const vector<Tc>& in_c);
    //operations
    template<typename Fc>
    void calculate_tcf(Fc fc, int i=0, int j=0);
    template<typename Fc>
    void calculate_tcf_all(Fc fc);
    
    void label_condition();
    
protected:
    template<typename Fc>
    void c_tcf_init(shared_ptr<vector<dipole_type> > v1, shared_ptr<vector<dipole_type> > v2,shared_ptr<vector<Tc> > c1, shared_ptr<vector<Tc> > c2,Fc fc,int nmax);
    template<typename Fc>
    void c_tcf(shared_ptr<vector<dipole_type> > v1, shared_ptr<vector<dipole_type> > v2,shared_ptr<vector<Tc> > c1, shared_ptr<vector<Tc> > c2,Fc fc,int nmax);
    //test
    bool allocate_condition(size_t i=0, size_t j=0) const;
};

//IR
inline bool IR::allocate_dt() const{
    return dt!=-1;
}
inline bool IR::allocate_dipole(size_t i) const{
    return dipole and i<dipole->size();
}
inline bool IR::allocate_vdipole(size_t i) const{
    return vdipole and i<vdipole->size();
}
inline bool IR::allocate_tcf(size_t i) const{
    return i<tcf.size();
}
inline bool IR::allocate_ir(size_t i) const{
    return i<ir.size();
}
inline bool IR::allocate_V() const{
    return V!=-1;
}
//IR2
inline bool IR2::allocate_dipole(size_t i, size_t j) const{
    return i<dipole.size() and dipole[i] and j<dipole[i]->size();
}
inline bool IR2::allocate_vdipole(size_t i, size_t j) const{
    return i<vdipole.size() and vdipole[i] and j<vdipole[i]->size();
}
//IR3
template<typename Tc>
IR3<Tc>::IR3(){
    //cout << " T: " << T << ", screen_t: " << screen_t << endl;
}
template<typename Tc>
void IR3<Tc>::set_condition(const vector<shared_ptr<vector<Tc> > >& in_c){
    condition = in_c;
}
template<typename Tc>
void IR3<Tc>::push_back_condition(shared_ptr<vector<Tc> > in_c){
    condition.push_back(in_c);
}
template<typename Tc>
void IR3<Tc>::push_back_condition(const vector<Tc>& in_c){
    condition.push_back(make_shared<vector<Tc>>(in_c));
}
template<typename Tc>
template<typename Fc>
void IR3<Tc>::calculate_tcf(Fc fc, int i, int j){
    if(!allocate_vdipole(i) or !allocate_vdipole(j))
        calculate_vdipole();
    test(allocate_vdipole(i) and allocate_vdipole(j),"calculate_tcf::no_vdipole_i_j");
    test(allocate_condition(i) and allocate_condition(j),"calculate_tcf::no_condition_i_j");
    os() << "Corelating " << setw(5) << i << " and " << setw(5) << j << '\r' << flush;
    if(allocate_tcf())
        c_tcf(vdipole.at(i),vdipole.at(j),condition.at(i),condition.at(j),fc,tcf_nmax);
    else
        c_tcf_init(vdipole.at(i),vdipole.at(j),condition.at(i),condition.at(j),fc,tcf_nmax);
}
template<typename Tc>
template<typename Fc>
void IR3<Tc>::calculate_tcf_all(Fc fc){
    for(int i=0;i<vdipole.size();i++){
        for(int j=0;j<vdipole.size();j++){
            calculate_tcf(fc,i,j);
        }
    }
}
template<typename Tc>
template<typename Fc>
void IR3<Tc>::c_tcf_init(shared_ptr<vector<dipole_type> > v1, shared_ptr<vector<dipole_type> > v2,shared_ptr<vector<Tc> > c1, shared_ptr<vector<Tc> > c2,Fc fc,int nmax){
    correlate3_function(v1->cbegin(),v1->cend(),v2->cbegin(),v2->cend(),back_inserter(tcf),c1->cbegin(),c2->cbegin(),fc,nmax,1);
}
template<typename Tc>
template<typename Fc>
void IR3<Tc>::c_tcf(shared_ptr<vector<dipole_type> > v1, shared_ptr<vector<dipole_type> > v2,shared_ptr<vector<Tc> > c1, shared_ptr<vector<Tc> > c2,Fc fc,int nmax){
    vector<double> tmp_tcf;
    correlate3_function(v1->cbegin(),v1->cend(),v2->cbegin(),v2->cend(),back_inserter(tmp_tcf),c1->cbegin(),c2->cbegin(),fc,nmax,1);
    auto tmp = tcf.begin();
    for( auto i : tmp_tcf)
        *(tmp++) += i;
}
template<typename Tc>
bool IR3<Tc>::allocate_condition(size_t i, size_t j) const{
    return i<condition.size() and condition[i] and j<condition[i]->size();
}

#endif
