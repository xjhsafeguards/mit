#ifndef ANWATER_H
#define ANWATER_H

#include "ccell.h"
class Water;

class Water{
    
    shared_ptr<Cell> cel;
    int indexO;
    vector<int> indexH;
    vector<double>  disH; // not used

public:
    static double OH_distance;
    
    Water() = default;
    Water(shared_ptr<Cell> in_cel);
    Water(shared_ptr<Cell> in_cel, int in_indexO);
    
    //get information
    int get_indexO() const;
    int nH() const;
    int get_indexH(int H) const;
    //const vector<int>& get_indexH() const;
    const typename Cell::position_type& posO() const;
    const typename Cell::position_type posO_incell() const;
    const typename Cell::position_type& posH(int H) const;
    const typename Cell::position_type& fposO() const;
    const typename Cell::position_type& fposH(int H) const;
    
    void write_index(ostream& os) const;
    
    //set information
    void set_indexO(int O);
    void push_back_H(int H);
    void clear_H();
    void set_indexH(const vector<int>& in_H);
    void set_cel(shared_ptr<Cell> in_cel);
    
    void set_water(); // assume "O" and "H" exist

protected:
    bool allocate_H(size_t i=0) const;
    void test(bool condition,string information) const;
};

//Water
inline int Water::get_indexO() const{
    return indexO;
}
inline int Water::nH() const{
    return indexH.size();
}
inline int Water::get_indexH(int H) const{
#ifndef N_TEST
    test(allocate_H(H),"get_indexH");
#endif
    return indexH.at(H);
}
//inline const vector<int>& Water::get_indexH() const{
//    return indexH;
//}
inline const typename Cell::position_type& Water::posO() const{
    return cel->apos("O",indexO);
}
inline const typename Cell::position_type Water::posO_incell() const{
    return cel->apos_incell("O",indexO);
}
inline const typename Cell::position_type& Water::posH(int H) const{
#ifndef N_TEST
    test(allocate_H(H),"posH");
#endif
    return cel->apos("H",indexH.at(H));
}
inline const typename Cell::position_type& Water::fposO() const{
    return cel->afpos("O",indexO);
}
inline const typename Cell::position_type& Water::fposH(int H) const{
#ifndef N_TEST
    test(allocate_H(H),"posH");
#endif
    return cel->afpos("H",indexH.at(H));
}
inline void Water::set_indexO(int O){
    indexO = O;
}
inline void Water::push_back_H(int H){
    indexH.push_back(H);
}
inline void Water::clear_H(){
    indexH.clear();
}
inline void Water::set_indexH(const vector<int>& in_H){
    indexH = in_H;
}
inline void Water::set_cel(shared_ptr<Cell> in_cel){
    cel = in_cel;
}
inline bool Water::allocate_H(size_t i) const{
    return i<indexH.size();
}

#endif
