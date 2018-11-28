#ifndef GDATA_H
#define GDATA_H

#include "gfun.h"
#include "gfile.h"
#include "gphy.h"
/*
 DATA TYPE:
 
 type:
    label_type;
 initialize:
 input:
    void set_label(is);
    //void check_label(is);
 output:
    label_type label();
 tool:
    void log();
    bool allocate_label();
    void test();
    void test_init();
    void further_throw(runtime_error err, string detail) const;
 */

/*
 POSITION_DATA TYPE:
 
 type:
    position_type;
 initialize:
    Position_data(string in_name, int in_na=-1);
 input:
    void set_na();
    void set_name();
 output:
    string get_name()
    int num();
    int get_na();
    bool isfractional();
    position_type [];
    iterator begin();
    iterator end();
    citerator cbegin() const;
    citerator cend() const;
 tool:
    bool allocate_{na(),name(),position()}
    void unset_position();
    void position_push_back(is,unit);
    void setposition(is,unit);
 
 */

/*
 Dictionary_DATA TYPE:
 
 type:
    key_value_type;
    mapped_value_type;
 initialize:
 input:
 output:
 
 tool:
    void dic_insert(is);
    sstream stream();
    bool {allocate,skip}_check(k);

 */
class Basic_Class;
class Data;
class Position_data;
class Dictionary_data;

class Basic_Class{
    
protected:
    // log information
    virtual void log(string information) const;
    virtual ostream& log() const;
    virtual ostream& os() const;
    // test
    virtual void test(bool condition, string error = "") const;
};

class Data: public Basic_Class{
    
public:
    typedef string label_type;
    
protected:
    shared_ptr<vector<label_type> >  labels =   make_shared<vector<label_type> >(); // label informations
    
public:
    //get information
    inline label_type label(int i=0) const;
    inline vector<label_type>&  get_labels() const;
    //set information
    inline void reset_label();
    void set_label(istream &is);
    void set_label(vector<label_type> in_label);
    virtual void set_label_file(Inputfile &inf);
    virtual void check_label(istream &is);
    
protected:
    // log information
    //virtual void log(string information) const;
    //virtual ostream& log() const;
    //virtual ostream& os() const;
    // test
    inline bool allocate_label(size_t i=0) const;
    virtual void test(bool condition, string error = "") const;
    virtual void test_init(bool initcondition, string name) const;
    void further_throw(runtime_error err, string detail) const;
};

class Position_data: public Data{
 
public:
    typedef Vector3<double> position_type;

    
protected:
    string          name    =   "default";      // name of the position
    int             na      =   -1;           // numbers of the position
    bool            fractional = false;       // whether fractional
    //shared_ptr<vector<string> >  labels =   make_shared<vector<string> >(); // label informations
    
    shared_ptr<vector<position_type>>   position   = make_shared<vector<Vector3<double> > >();

public:
    typedef decltype(position->begin())     iterator;
    typedef decltype(position->cbegin())    citerator;
    Position_data() = default;
    Position_data(string in_name, int in_na=-1);
    //operation
    virtual void value_duplicate();
    //get information
    string get_name() const;
    int num() const;
    int get_na() const;
    bool isfractional() const;
    position_type& operator[](int i);
    const position_type& operator[](int i) const;
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
    
private:
    position_type buffer;   // used to read data
};

/*
class Dictionary_data: public Data{
    
public:
    typedef string key_value_type;
    typedef string mapped_value_type;
    
protected:
    key_value_type skiper = "!";
    shared_ptr<map<key_value_type,mapped_value_type> >  dic =   make_shared<map<key_value_type,mapped_value_type> >();
    
public:
    
protected:
    //test
    bool allocate_check(key_value_type k) const;
    bool skip_check(key_value_type k) const;
    //set information
    void dic_insert(istream& is);
    void dic_insert(istream& is, key_value_type k);
    //read information
    stringstream& stream(key_value_type k) const;
    
private:
    key_value_type    key_buffer;  // used to read data
    mapped_value_type map_buffer;  // used to read data
    mutable stringstream buffer;
};
 */

//basic_class
inline void Basic_Class::log(string information) const {
    ofs_log << " " << information << endl;
}

inline ostream& Basic_Class::log() const {
    return ofs_log;
}
inline ostream& Basic_Class::os() const {
    return cout;
}
inline void Basic_Class::test(bool condition, string error) const{
    if(!condition)
        throw(runtime_error("Fail_"+error));
}
//Data
inline typename Data::label_type Data::label(int i) const{
#ifndef N_TEST
    test(allocate_label(i), "label()::out_of_range");
#endif
    return (*labels)[i];
}
inline vector<typename Data::label_type>&  Data::get_labels() const{
    return *labels;
}
inline void Data::reset_label(){
    labels =   make_shared<vector<label_type> >();
}
inline bool Data::allocate_label(size_t i) const {
    return i<labels->size();
}

inline void Data::test(bool condition, string error) const {
    if(!condition)
        throw(runtime_error("DATA::Fail_"+error));
}
inline void Data::test_init(bool initcondition, string name) const {
    if(initcondition)
        throw(runtime_error("DATA::Wrong_initialze:" + name));
}
inline void Data::further_throw(runtime_error err, string detail) const{
    throw(runtime_error(static_cast<string>(err.what())+ detail));
}
//Position_data
inline string Position_data::get_name() const {
#ifndef N_TEST
    test_init(!allocate_name(),name+":na");
#endif
    return name;
}
inline int Position_data::num() const {
    return position->size();
}
inline int Position_data::get_na() const {
#ifndef N_TEST
    test_init(!allocate_na(),name+":na");
#endif
    return na;
}
inline bool Position_data::isfractional() const {
    return fractional;
}
inline typename Position_data::position_type& Position_data::operator[](int i){
#ifndef N_TEST
    test(allocate_position(i), name+"::operator[]::out_of_range");
#endif
    return (*position)[i];
}
inline const typename Position_data::position_type& Position_data::operator[](int i) const {
#ifndef N_TEST
    test(allocate_position(i), name+"::operator[]::out_of_range");
#endif
    return (*position)[i];
}
inline typename Position_data::iterator Position_data::begin(){
    return position->begin();
}
inline typename Position_data::iterator Position_data::end(){
    return position->end();
}
inline typename Position_data::citerator Position_data::cbegin() const{
    return position->cbegin();
}
inline typename Position_data::citerator Position_data::cend() const{
    return position->cend();
}
inline void Position_data::set_na(int in_na){
    na = in_na;
}
inline void Position_data::set_name(string in_name){
    name = in_name;
}
inline void Position_data::set_fractional(bool fra){
    fractional = fra;
}
inline bool Position_data::allocate_na() const {
    return (na > -1);
}
inline bool Position_data::allocate_name() const {
    return (name != "default");
}
inline bool Position_data::allocate_position(size_t i) const {
    return i<position->size();
}
inline void Position_data::unset_position(){
    position->clear();
}
inline void Position_data::position_push_back(istream& is,double unit_conv){
    is >> buffer;
    is.ignore(500,'\n');
    position->push_back(buffer*unit_conv);
}
inline void Position_data::skip_position(istream& is) const{
#ifndef N_TEST
    test(is.good(),"Skip_position");
#endif
    for(int i=0;i<na;i++)
        is.ignore(5000,'\n');
}
/*
//Dictionary_data
inline bool Dictionary_data::allocate_check(key_value_type k) const{
    return dic->find(k) != dic->end();
}
inline bool Dictionary_data::skip_check(key_value_type k) const{
    return k.find(skiper) == string::npos;
}
inline void Dictionary_data::dic_insert(istream& is){
    is >> key_buffer;
    getline(is,map_buffer);
    if(allocate_check(key_buffer))
        log(" Discard duplated input:" + key_buffer);
    else if(!skip_check(key_buffer))
        dic_insert(is,key_buffer);
}
inline void Dictionary_data::dic_insert(istream& is, key_value_type k){
    getline(is,map_buffer);
    if(!isspace(map_buffer))
        dic->insert(make_pair(key_buffer,map_buffer));
    else
        log(" Discard empty input:" + k);
}
inline stringstream& Dictionary_data::stream(key_value_type k) const{
    buffer.str("");
    test(allocate_check(k),"DICTINOARY::STREAM::check_key:" + k);
    buffer << dic->at(k);
    return buffer;
}*/

#endif
