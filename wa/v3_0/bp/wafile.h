#ifndef WAFILE_H
#define WAFILE_H

class Wafile;
class Wainfile;
class Waoutfile;
class Waiofile;

//#include "gfun.h"

class Wafile{
    
protected:
    string filename = "noname";
    string filetype = "notype";
public:
    //constructor
    Wafile() = default;
    Wafile(string name);
    Wafile(string name,string type);
    
    //get information
    string name() const;
    string name_only() const;
    string type() const;
    //set files
    void set_name(string name);
    void set_type();
    void set_type(string type);
    virtual void open() {}
    
protected:
    //test
    bool allocate_name() const;
    bool allocate_type() const;
    void test(bool condition, string error = "") const;
    void set_suffix_filetype();

    
};
