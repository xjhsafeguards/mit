#ifndef WAHBOND_H
#define WAHBOND_H

#include "anhbond.h"
#include "wawater.h"

class Wahbonds: public Waters{

protected:
    vector<shared_ptr<Hbond> >  hbonds;
    //Hbonds information for those water with numH != 2
    double avg_aion=0;
    double avg_dion=0;
    double avg_count_ion=0;
    
    //Hbonds information for those water with numH == 2
    double avg_awater=0;
    double avg_dwater=0;
    double avg_count_water=0;
    
public:
    static double OO_distance;//to find H_bond(SI)
    static double HOO_angle;//to find H_bond(degree)
    static void Routine();
    
    Wahbonds() = default;
    Wahbonds(Cell& in_cel, bool ifsave=0);
    //get information
    int nhbond(int i) const;
    int naccept(int i) const;
    int ndonate(int i) const;
    const vector<int>& get_index_accept(int i) const;
    const vector<int>& get_index_donate(int i) const;
    const vector<double>&  get_donate_angle(int i) const;
    int get_index_accept(int i,int ind) const;
    int get_index_donate(int i,int ind) const;
    double get_donate_angle(int i,int ind) const;
    bool check_accept(int j,int i) const;
    bool check_donate(int j,int i) const;
    bool check_bond(int j,int i) const;
    
    static void write_average_hbonds_header(ostream& os);
    void write_average_hbonds(ostream& os) const;
    static void write_each_hbonds_header(ostream& os);
    void write_each_hbonds(ostream& os) const;
    static void write_ions_hbonds_header(ostream& os);
    void write_ions_hbonds(ostream& os) const;
    
    //set information
    void init_hbonds();
    void set_hbonds();
    void unset_hbonds();
    
private:
    //get
    shared_ptr<Hbond> get_hbonds(int i) const;
    //set
    void donate(int O1,int O2);
    //test
    bool allocate_hbonds(size_t i=0) const;
    void test(bool condition,string information) const;
};

class WahbondsT{
    
    int count=0;
    int countion=0;
    map<int,int> count_iona;
    map<int,int> count_iond;
    map<int,int> count_water;
    int totaion=0;
    int totdion=0;
    int tothwater=0;
    
    //int count_iona[10]={0};
    //int count_iond[10]={0};
    //int count_water[10]={0};
    
public:
    
    static void Routine();
    WahbondsT() = default;
    //get
    void write_water(ostream& os);
    void write_ion(ostream& os);
    
    //set
    void add(const Wahbonds& WF);
    void add_water(const Wahbonds& WF);
    void add_water(const Wahbonds& WF, double z1, double z2);
    
private:
    //set
    void set_totaion();
    void set_totdion();
    void set_tothwater();
};

inline int Wahbonds::nhbond(int i) const{
#ifndef N_TEST
    test(allocate_hbonds(i),"nhb");
#endif
    return hbonds[i]->n();
}
inline int Wahbonds::naccept(int i) const{
#ifndef N_TEST
    test(allocate_hbonds(i),"naccept");
#endif
    return hbonds[i]->naccept();
}
inline int Wahbonds::ndonate(int i) const{
#ifndef N_TEST
    test(allocate_hbonds(i),"ndonate");
#endif
    return hbonds[i]->ndonate();
}
inline const vector<int>& Wahbonds::get_index_accept(int i) const{
    return get_hbonds(i)->get_index_accept();
}
inline const vector<int>& Wahbonds::get_index_donate(int i) const{
    return get_hbonds(i)->get_index_donate();
}
inline const vector<double>& Wahbonds::get_donate_angle(int i) const{
    return get_hbonds(i)->get_donate_angle();
}
inline int Wahbonds::get_index_accept(int i,int ind) const{
    return get_hbonds(i)->get_index_accept(ind);
}
inline int Wahbonds::get_index_donate(int i,int ind) const{
    return get_hbonds(i)->get_index_donate(ind);
}
inline double Wahbonds::get_donate_angle(int i,int ind) const{
    return get_hbonds(i)->get_donate_angle(ind);
}
inline bool  Wahbonds::check_accept(int j,int i) const{
#ifndef N_TEST
    test(allocate_hbonds(j),"check_accept");
#endif
    return hbonds[j]->check_accept(i);
}
inline bool  Wahbonds::check_donate(int j,int i) const{
#ifndef N_TEST
    test(allocate_hbonds(j),"check_donate");
#endif
    return hbonds[j]->check_donate(i);
}
inline bool  Wahbonds::check_bond(int j,int i) const{
#ifndef N_TEST
    test(allocate_hbonds(j),"check_bond");
#endif
    return hbonds[j]->check_bond(i);
}
inline void Wahbonds::unset_hbonds(){
    hbonds.clear();
}
inline shared_ptr<Hbond> Wahbonds::get_hbonds(int i) const{
#ifndef N_TEST
    test(allocate_hbonds(i),"get_hbonds");
#endif
    return hbonds[i];
}
inline bool Wahbonds::allocate_hbonds(size_t i) const{
    return i<hbonds.size() and hbonds[i];
}
#endif
