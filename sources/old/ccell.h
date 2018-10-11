#ifndef CCELL_H
#define CCELL_H

#include "gheader.h"
#include "input.h"
#include "catom.h"
#include "ccellp.h"
#include "cwannier.h"

/*
 CELL DATA
 type:
    vector_type;
    position_type;
 initialize:
    Cell(int ntype, bool in_cellp);
    Cell(int ntype, bool in_cellp, int nwan);
 */

class Cell;
Inputfile& operator>>(Inputfile&, Cell&);
Outputfile& operator<<(Outputfile&, const Cell&);
Outputfile& operator<<(Outputfile&,Cell&);

class Cell: public Data{
    friend Inputfile& operator>>(Inputfile&, Cell&);
    friend Outputfile& operator<<(Outputfile&, const Cell&);
    friend Outputfile& operator<<(Outputfile&,Cell&);
public:
    typedef typename Cellp::vector_type         vector_type;
    typedef typename Cellp::value_type          cellp_value_type;
    typedef typename Position_data::position_type position_type;
    typedef typename Position_data::iterator    piterator;
    typedef typename Position_data::citerator   pciterator;
    
protected:
    
    shared_ptr<Cellp>           cellp;
    shared_ptr<vector<Atom> >   atoms;
    shared_ptr<Wannier>         wannier;
    //shared_ptr<vector<string> >  labels =   make_shared<vector<string> >(); // label informations
    // 0 snapshot 1 time;
    shared_ptr<vector<Atom> >   fatoms;
    shared_ptr<Wannier>         fwannier;
    
    //int         celltype = 0;   //
    Inputfile   geo_file;
    Inputfile   cel_file;
    Inputfile   wan_file;

    static unsigned start;
    static unsigned step;
    static unsigned end;
    unsigned count  =0;
    unsigned num    =0;
    double volume_sum = 0;
    
    static map<string,int> index_dic;
    
public:
    Cell() = default;
    Cell(Input input);
    //operation
    void init();
    bool next();
    bool inrange() const;
    unsigned get_num() const;
    unsigned get_count() const;
    //get information
    piterator abegin(int i);
    piterator aend(int i);
    pciterator acbegin(int i) const;
    pciterator acend(int i) const;
    piterator abegin(string name);
    piterator aend(string name);
    pciterator acbegin(string name) const;
    pciterator acend(string name) const;
    piterator wbegin();
    piterator wend();
    pciterator wcbegin() const;
    pciterator wcend() const;
    const position_type& apos(int i,int j) const;
    const position_type& apos(string name,int j) const;
    const position_type& wpos(int i) const;
    const position_type apos_incell(int i,int j) const;
    const position_type apos_incell(string name,int j) const;
    const position_type wpos_incell(int i) const;
    const position_type& afpos(int i,int j) const;
    const position_type& afpos(string name,int j) const;
    const position_type& wfpos(int i) const;
    const cellp_value_type& box() const;
    int anum() const;// total atom numbers allocated
    int anum(int i) const;// ith atom numbers
    int anum(string name) const;
    int ana() const;
    int ana(int i) const;
    int ana(string name) const;
    int aindex(string name) const;
    int wnum() const;
    int wna() const;
    int snapshot() const;
    double time() const;
    double volume() const;
    double avolume() const;
    shared_ptr<Cell> save() const; //return a unchanged cel while cel reads
    

    //operation
    static position_type Cartesian(const cellp_value_type& cp,const position_type& a);
    static position_type In_cell(const position_type& a);
    static position_type Shortest_vector(const position_type& a,const position_type& b);
    static double Distance(const cellp_value_type& cp, const position_type& a, const position_type& b);
    
    double distance(int i1,int j1, int i2, int j2) const;
    double distance(string n1,int j1,string n2, int j2) const;
    double distance(int i1,int j1, int w2) const;
    double distance(string n1,int j1, int w2) const;
    double angle(int i1,int j1, int i2, int j2, int i3, int j3) const;
    double angle(string n1,int j1,string n2, int j2,string n3, int j3) const;
    vector_type cvector(int i1,int j1, int i2, int j2) const;
    vector_type cvector(string n1,int j1,string n2, int j2) const;
    vector_type cvector(int i1,int j1, int w2) const;
    vector_type cvector(string n1,int j1, int w2) const;
    //set information
    void set();
    void skip();
    
    // IO
    void read_cif(istream&);
    
    void write_in(ostream&) const;
    void write_POSCAR(ostream&) const;
    
protected:
    //set
    void set_atoms_nfractional();
    void set_wannier_nfractional();
    void set_fatoms_fractional();
    void set_fwannier_fractional();
    void set_up_fatoms(); // set up and fatoms after allocate atoms
    void set_geo();
    void set_wan();
    void set_cel();
    void skip_geo();
    void skip_wan();
    void skip_cel();
    void set_pos_friction(Position_data& pos);
    void set_pos_cartesian(Position_data& pos);
    void check_labels(vector<string>);
    
    void set_read_cif(istream&); // allocate cellp atoms and set na and name for atoms
    
    //operation for fractional position
    position_type cartesian(const position_type& a) const;
    position_type in_cell(const position_type& a) const;
    position_type shortest_vector(const position_type& a,const position_type& b)const;// return from a to b;
    //test
    inline bool allocate_atoms(size_t i=0) const;
    inline bool allocate_fatoms(size_t i=0) const;
    inline bool allocate_wannier() const;
    inline bool allocate_fwannier() const;
    inline bool allocate_cellp() const;
};

inline bool Cell::inrange() const{
    return num+step<=end;
}
inline unsigned Cell::get_num() const{
    return num;
}
inline unsigned Cell::get_count() const{
    return count;
}
inline typename Cell::piterator Cell::abegin(int i){
#ifndef N_TEST
    test(allocate_atoms(i),"Cell::abegin::out_of_range");
#endif
    return atoms->at(i).begin();
}
inline typename Cell::piterator Cell::aend(int i){
#ifndef N_TEST
    test(allocate_atoms(i),"Cell::abegin::out_of_range");
#endif
    return atoms->at(i).end();
}
inline typename Cell::pciterator Cell::acbegin(int i) const{
#ifndef N_TEST
    test(allocate_atoms(i),"Cell::abegin::out_of_range");
#endif
    return atoms->at(i).cbegin();
}
inline typename Cell::pciterator Cell::acend(int i) const{
#ifndef N_TEST
    test(allocate_atoms(i),"Cell::abegin::out_of_range");
#endif
    return atoms->at(i).cend();
}
inline typename Cell::piterator Cell::abegin(string name){
    return atoms->at(aindex(name)).begin();
}
inline typename Cell::piterator Cell::aend(string name){
    return atoms->at(aindex(name)).end();
}
inline typename Cell::pciterator Cell::acbegin(string name) const{
    return atoms->at(aindex(name)).cbegin();
}
inline typename Cell::pciterator Cell::acend(string name) const{
    return atoms->at(aindex(name)).cend();
}
inline typename Cell::piterator Cell::wbegin(){
#ifndef N_TEST
    test(allocate_wannier(),"Cell::wpos::no_wannier");
#endif
    return wannier->begin();
}
inline typename Cell::piterator Cell::wend(){
#ifndef N_TEST
    test(allocate_wannier(),"Cell::wpos::no_wannier");
#endif
    return wannier->end();
}
inline typename Cell::pciterator Cell::wcbegin() const{
#ifndef N_TEST
    test(allocate_wannier(),"Cell::wpos::no_wannier");
#endif
    return wannier->cbegin();
}
inline typename Cell::pciterator Cell::wcend() const{
#ifndef N_TEST
    test(allocate_wannier(),"Cell::wpos::no_wannier");
#endif
    return wannier->cend();
}
inline const typename Cell::position_type& Cell::apos(int i,int j) const{
#ifndef N_TEST
    test(allocate_atoms(i),"Cell::apos::out_of_range");
#endif
    return atoms->at(i)[j];
}
inline const typename Cell::position_type& Cell::apos(string name,int j) const{
    return atoms->at(aindex(name))[j];
}
inline const typename Cell::position_type& Cell::wpos(int i) const{
#ifndef N_TEST
    test(allocate_wannier(),"Cell::wpos::no_wannier");
#endif
    return (*wannier)[i];
}
inline const typename Cell::position_type Cell::apos_incell(int i,int j) const{
#ifndef N_TEST
    test(allocate_atoms(i),"Cell::apos::out_of_range");
#endif
    return cartesian(in_cell(fatoms->at(i)[j]));
}
inline const typename Cell::position_type Cell::apos_incell(string name,int j) const{
    return cartesian(in_cell(fatoms->at(aindex(name))[j]));
}
inline const typename Cell::position_type Cell::wpos_incell(int i) const{
#ifndef N_TEST
    test(allocate_wannier(),"Cell::wpos::no_wannier");
#endif
    return cartesian(in_cell((*fwannier)[i]));
}
inline const typename Cell::position_type& Cell::afpos(int i,int j) const{
#ifndef N_TEST
    test(allocate_atoms(i),"Cell::afpos::out_of_range");
#endif
    return fatoms->at(i)[j];
}
inline const typename Cell::position_type& Cell::afpos(string name,int j) const{
    return fatoms->at(aindex(name))[j];
}
inline const typename Cell::position_type& Cell::wfpos(int i) const{
#ifndef N_TEST
    test(allocate_wannier(),"Cell::wfpos::no_wannier");
#endif
    return (*fwannier)[i];
}
inline const typename Cell::cellp_value_type& Cell::box() const{
#ifndef N_TEST
    test(allocate_cellp(),"Cell::box::no_cellp");
#endif
    return (*cellp)();
}
inline int Cell::snapshot() const{
    if(allocate_label(0))
        return stoi(label(0));
    return 0;
}
inline double Cell::time() const{
    if(allocate_label(1))
        return stod(label(1));
    return 0;
}
inline double Cell::volume() const{
#ifndef N_TEST
    test(allocate_cellp(),"Cell::volume::no_cellp");
#endif
    return (*cellp).volume();
}
inline void Cell::set(){
    reset_label();
    set_cel();
    set_geo();
    set_wan();
    volume_sum += volume();
    ++num;
    os() << setw(16) << "read snapshot: " << setw(10) << num << setw(3) << num*100/end << "%"<< '\r' << flush;
    count++;
}
inline void Cell::skip(){
    skip_cel();
    skip_geo();
    skip_wan();
    ++num;
    os() << setw(16) << "skip snapshot: " << setw(10) << num << setw(3) << num*100/end << "%"<< '\r' << flush;
}
inline typename Cell::position_type Cell::cartesian(const position_type& a) const{
    return a*box();
}
inline typename Cell::position_type Cell::in_cell(const position_type& a) const{
    position_type b(a);
    b.normal_BC(position_type(1.0,1.0,1.0));
    return b;
}
inline typename Cell::position_type Cell::shortest_vector(const position_type& a,const position_type& b) const{
    position_type c = b - a;
    c = in_cell(c+position_type(0.5,0.5,0.5)) - position_type(0.5,0.5,0.5);
    return c;
}
inline bool Cell::allocate_atoms(size_t i) const{
    return i<atoms->size();
}
inline bool Cell::allocate_fatoms(size_t i) const{
    return i<fatoms->size();
}
inline bool Cell::allocate_wannier() const{
    return static_cast<bool>(wannier);
}
inline bool Cell::allocate_fwannier() const{
    return static_cast<bool>(fwannier);
}
inline bool Cell::allocate_cellp() const{
    return static_cast<bool>(cellp);
}
//Cell(int ntype);
//Cell(int ntype,int nwan);

#endif
