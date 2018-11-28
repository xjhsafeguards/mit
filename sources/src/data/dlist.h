#ifndef DLIST_H
#define DLIST_H

#include "gheader.h"
#include "deigen.h"

//Dlist a value like class
class Dlist{
    int              na = -1;
//    double           mass = -1;
//    double           charge = -1;
//    string           name = "default";
    vector<vec_type> positions;
    bool             fractional = false;
public:
    typedef decltype(positions.begin())     iterator;
    typedef decltype(positions.cbegin())    citerator;

    Dlist() = default;
    Dlist(int in_na=-1);

    string get_name() const;
    int num() const;
    int get_na() const;
    bool isfractional() const;
    vec_type& operator[](int i);
    const vec_type& operator[](int i) const;
    iterator begin();
    iterator end();
    citerator cbegin() const;
    citerator cend() const;

    //set information
    void set_na(int in_na);
    void set_name(string in_name);
    void set_fractional(bool fra);

protected:
    //test
    bool allocate_na() const;
    bool allocate_name() const;
    bool allocate_position(size_t i=0) const;
    void unset_position();
    void position_push_back(istream& is,double unit_conv=1);
    void set_position(istream& is,double unit_conv=1);
    void skip_position(istream& is) const;
};

#endif
