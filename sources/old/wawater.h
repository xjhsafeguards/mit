#ifndef WAWATER_H
#define WAWATER_H

#include "anwater.h"
class Waters;
Outputfile& operator<<(Outputfile&,Waters&);

class Waters{
    friend Outputfile& operator<<(Outputfile&,Waters&);
    
protected:
    shared_ptr<Cell>             cel;
    vector<shared_ptr<Water> >   waters;        // include ions
    vector<int>                  ions;
    
public:
    Waters(Cell& in_cel,bool ifsave=0);
    //get information
    const Cell& get_cell() const;
    const Water& get_water(int i) const;
    int get_indexO(int i)const;
    int get_indexH(int i,int j) const;
    int get_nH(int i) const;
    const typename Cell::position_type& posO(int i) const;
    const typename Cell::position_type posO_incell(int i) const;
    const typename Cell::position_type& posH(int i,int H) const;
    const typename Cell::position_type& fposO(int i) const;
    const typename Cell::position_type& fposH(int i,int H) const;
    int nwater() const;
    int get_ion(int i) const;
    int nion() const;
    double time() const;
    int snapshot() const;
    
    void write_waters(ostream& os);
    
    //set information
    void set_waters();
    void unset_waters();
    void unset_cel();
    
protected:
    bool allocate_waters(size_t i=0) const;
    bool allocate_ions(size_t i=0) const;
    void test(bool condition,string information) const;
    void log(string information) const;
    ostream& log() const;
};

//Waters
inline const Cell& Waters::get_cell() const{
    //test(cel,"get_cell::no_cel");
    return *cel;
}
inline const Water& Waters::get_water(int i) const{
#ifndef N_TEST
    test(allocate_waters(i),"get_water");
#endif
    return *(waters.at(i));
}
inline int Waters::get_indexO(int i)const{
    return get_water(i).get_indexO();
}
inline int Waters::get_indexH(int i,int j) const{
    return get_water(i).get_indexH(j);
}
inline int Waters::get_nH(int i) const{
    return get_water(i).nH();
}
inline const typename Cell::position_type& Waters::posO(int i) const{
    return get_water(i).posO();
}
inline const typename Cell::position_type Waters::posO_incell(int i) const{
    return cel->apos_incell("O",get_water(i).get_indexO());
}
inline const typename Cell::position_type& Waters::posH(int i,int H) const{
    return get_water(i).posH(H);
}
inline const typename Cell::position_type& Waters::fposO(int i) const{
    return get_water(i).fposO();
}
inline const typename Cell::position_type& Waters::fposH(int i,int H) const{
    return get_water(i).fposH(H);
}
inline int Waters::nwater() const{
    return waters.size();
}
inline int Waters::get_ion(int i) const{
#ifndef N_TEST
    test(allocate_ions(i),"get_ion");
#endif
    return ions[i];
}
inline int Waters::nion() const{
    return ions.size();
}
inline double Waters::time() const{
    return cel->time();
}
inline int Waters::snapshot() const{
    return cel->snapshot();
}
inline void Waters::unset_waters(){
    waters.clear();
    ions.clear();
}
inline void Waters::unset_cel(){
    cel.reset();
}
inline bool Waters::allocate_waters(size_t i) const{
    return i<waters.size() and waters[i];
}
inline bool Waters::allocate_ions(size_t i) const{
    return i<ions.size();
}
inline void Waters::log(string information) const{
    GLOG->stream() << " " << information << endl;
}
inline ostream& Waters::log() const{
    return GLOG->stream() << " ";
}
#endif
