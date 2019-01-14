#include "dlist.h"

Dlist::Dlist(int in_na=-1):na(in_na){

}

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

//test
bool allocate_na() const;
bool allocate_name() const;
bool allocate_position(size_t i=0) const;
void unset_position();
void position_push_back(istream& is,double unit_conv=1);
void set_position(istream& is,double unit_conv=1);
void skip_position(istream& is) const;
