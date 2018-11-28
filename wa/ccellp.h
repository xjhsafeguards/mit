#ifndef CCELLP_H
#define CCELLP_H

class Cellp;
#include "gheader.h"

/*
 CELLP DATA
 
 type:
    value_type
    vector_type
 initialize:
    Cellp(value_type in_mat);
 input:
    void read_cel(istream& is);
    void read_Input();
 output:
    vector_typ [];
    value_type operator()
 
 
 */
Inputfile& operator>>(Inputfile &,Cellp &);
Outputfile& operator<<(Outputfile&, const Cellp&);

//only work on mat
istream& operator>>(istream&,Cellp&);
ostream& operator<<(ostream&,const Cellp&);

class Cellp: public Data {
    friend Inputfile& operator>>(Inputfile &,Cellp &);
    friend Outputfile& operator<<(Outputfile&,const Cellp&);
    friend istream& operator>>(istream&,Cellp&);
    friend ostream& operator<<(ostream&,const Cellp&);
    friend class Cell;
public:
    typedef Matrix3<double> value_type;
    typedef Vector3<double> vector_type;
    
protected:
    shared_ptr<value_type>    mat;
    
public:
    Cellp() = default;
    Cellp(value_type in_mat);
    
    //get information
    inline vector_type& operator[](int i);
    inline const vector_type& operator[](int i) const;
    inline value_type& operator()();
    inline const value_type& operator()() const;
    inline double volume() const;
    shared_ptr<Cellp> save() const;
    void get_parameter(double&,double&,double&,double&,double&,double&) const;
    
    void write_mat(ostream&, double unitconv=1) const;
    inline void write_in(ostream&) const;
    inline void write_POSCAR(ostream&) const;
    
    //set information
    void set_mat(istream& is,double unitconv=1);
    void skip(Inputfile& Inf) const;
    void set_parameter(double,double,double,double,double,double,double unitconv=1);
    
    void read_cif(istream& is);
    inline void read_cel(istream& is);
    inline void skip_cel(istream& is) const;
    void read_Input();
    
protected:
    //test
    inline bool allocate_mat() const;
    //set information

};

typename Cellp::vector_type& Cellp::operator[](int i){
#ifndef N_TEST
    test(allocate_mat(),"Cellp::operator[]::no_mat");
#endif
    return (*mat)[i];
}
const typename Cellp::vector_type& Cellp::operator[](int i) const {
#ifndef N_TEST
    test(allocate_mat(),"Cellp::operator[]::no_mat");
#endif
    return (*mat)[i];
}
typename Cellp::value_type& Cellp::operator()(){
#ifndef N_TEST
    test(allocate_mat(),"Cellp::operator[]::no_mat");
#endif
    return *mat;
}
const typename Cellp::value_type& Cellp::operator()() const {
#ifndef N_TEST
    test(allocate_mat(),"Cellp::operator[]::no_mat");
#endif
    return *mat;
}
inline double Cellp::volume() const{
#ifndef N_TEST
    test(allocate_mat(),"Cellp::volume::no_mat");
#endif
    return (*mat)[0].cross((*mat)[1])*(*mat)[2];
}
inline void Cellp::write_in(ostream& os) const{
    write_mat(os);
}
inline void Cellp::write_POSCAR(ostream& os) const{
    write_mat(os);
}
inline void Cellp::read_cel(istream& is){
    set_label(is);
    set_mat(is,Unit::Bohr2A);
    is.ignore(500,'\n');
}
inline void Cellp::skip_cel(istream& is) const{
    for(int i=0; i!=4;i++)
        is.ignore(500,'\n');
}
inline bool  Cellp::allocate_mat() const{
    return static_cast<bool>(mat);
}
#endif
