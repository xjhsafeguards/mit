#ifndef ANHBOND_H
#define ANHBOND_H

#include "ccell.h"
class Hbond;

class Hbond{
  
    vector<int> index_accept;
    vector<int> index_donate;
    vector<double> donate_angle;
    
public:
    Hbond() = default;
    
    //get information
    int n() const;
    int naccept() const;
    int ndonate() const;
    int ndonate_angle() const;
    const vector<int>& get_index_accept() const;
    const vector<int>& get_index_donate() const;
    const vector<double>&  get_donate_angle() const;
    int get_index_accept(int i) const;
    int get_index_donate(int i) const;
    double get_donate_angle(int i) const;
    bool check_accept(int i) const;
    bool check_donate(int i) const;
    bool check_bond(int i) const;
    
    static void write_detail_header(ostream& os);
    void write_detail(ostream& os) const;
    
    //set information
    void push_back_accept(int H);
    void push_back_donate(int H);
    void push_back_angle(double ang);
    void clear();
    
private:
    bool allocate_accept(size_t i=0) const;
    bool allocate_donate(size_t i=0) const;
    bool allocate_angle(size_t i=0) const;
    void test(bool condition,string information) const;
};

inline int Hbond::n() const{
    return index_accept.size() + index_donate.size();
}
inline int Hbond::naccept() const{
    return index_accept.size();
}
inline int Hbond::ndonate() const{
    return index_donate.size();
}
inline int Hbond::ndonate_angle() const{
    return donate_angle.size();
}
inline const vector<int>& Hbond::get_index_accept() const{
    return index_accept;
}
inline const vector<int>& Hbond::get_index_donate() const{
    return index_donate;
}
inline const vector<double>&  Hbond::get_donate_angle() const{
    return donate_angle;
}
inline int Hbond::get_index_accept(int i) const{
#ifndef N_TEST
    test(allocate_accept(i),"get_index_accept");
#endif
    return index_accept.at(i);
}
inline int Hbond::get_index_donate(int i) const{
#ifndef N_TEST
    test(allocate_donate(i),"get_index_accept");
#endif
    return index_donate.at(i);
}
inline double Hbond::get_donate_angle(int i) const{
#ifndef N_TEST
    test(allocate_angle(i),"get_index_accept");
#endif
    return donate_angle.at(i);
}
inline void Hbond::push_back_accept(int H){
    index_accept.push_back(H);
}
inline void Hbond::push_back_donate(int H){
    index_donate.push_back(H);
}
inline void Hbond::push_back_angle(double ang){
    donate_angle.push_back(ang);
}
inline void Hbond::clear(){
    index_accept.clear();
    index_donate.clear();
    donate_angle.clear();
}
inline bool Hbond::allocate_accept(size_t i) const{
    return i<index_accept.size();
}
inline bool Hbond::allocate_donate(size_t i) const{
    return i<index_donate.size();
}
inline bool Hbond::allocate_angle(size_t i) const{
    return i<donate_angle.size();
}
#endif
