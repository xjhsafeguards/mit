#ifndef GFILE_H
#define GFILE_H

class File;
class Inputfile;
class Outputfile;

#include "gfun.h"

istream& operator>>(istream&, File&); // read in file names as a string
ostream& operator<<(ostream&, File&); // print name

class File
{
    
    friend istream& operator>>(istream&, File&);
    friend ostream& operator<<(ostream&, File&);
protected:
    
    string filename = "noname";
    string filetype = "notype";
public:
    //constructor
    File() = default;
    File(string name);
    File(string name,string type);
    
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
    //running log
    void log(string information) const;
    ostream& os(ostream& fos = cout) const;
    //test
    bool allocate_name() const;
    bool allocate_type() const;
    void test(bool condition, string error = "") const;
    void set_suffix_filetype();

};

class Inputfile: public File
{

    shared_ptr<ifstream> ifs;
    
public:
    
    Inputfile() = default;
    Inputfile(string );
    Inputfile(string,string );
    ~Inputfile() ;
    //set files
    void open();
    void open(string name);
    void seek_top(bool show_log = true);

    //tests
    operator bool() const;
    bool allocate() const;

    //read operations
    template <typename T> const
    Inputfile& operator>>(T& holder) const;
    ifstream& stream() const;
    
    
protected:
    bool allocate_ifs() const;
    bool good() const;
};

class Outputfile: public File
{
    
    shared_ptr<ofstream> ofs;
    
public:
    
    Outputfile() = default;
    Outputfile(string );
    Outputfile(string,string );
    ~Outputfile() ;
    //set files
    void open();
    void open(string name);
    void setprecision(size_t i=10);
    
    //tests
    operator bool() const;
    bool allocate() const;
    
    //read operations
    template <typename T> const
    Outputfile& operator<<(T& holder) const;
    ofstream& stream() const;
    
    
protected:
    bool allocate_ofs() const;
    bool good() const;
};
//File
inline string File::name() const {
    return filename;
}
inline string File::type() const {
    return filetype;
}
inline void File::set_name(string name){
    filename = name;
}
inline void File::set_type() {
    set_suffix_filetype();
}
inline void File::set_type(string type) {
    filetype = type;
}
inline void File::log(string information) const {
    GLOG->stream() << " " << information << endl;
}
inline ostream& File::os(ostream& fos) const{
    return fos;
}
inline bool File::allocate_name() const {
    return ((filename == "noname") ? false : true);
}
inline bool File::allocate_type() const {
    return ((filetype == "notype") ? false : true);
}
inline void File::test(bool condition, string error) const {
    if(!condition)
        throw(runtime_error("FILE:"+filename+"::Fail_"+error));
}

//Inputfile
inline Inputfile::operator bool() const {
    return good();
}
inline bool Inputfile::allocate() const {
    return allocate_ifs();
}
template <typename T> const
Inputfile& Inputfile::operator>>(T& holder) const {
    test(good(),"operator>>");
    *ifs >> holder;
    return *this;
}
inline ifstream& Inputfile::stream() const {
    test(good(),"stream");
    return *ifs;
}
inline bool Inputfile::allocate_ifs() const {
    return static_cast<bool>(ifs);
}
inline bool Inputfile::good() const {
    return (allocate_ifs() ? ifs->good() : false);
}

//Outputfile
inline void Outputfile::setprecision(size_t i){
    (*ofs).precision(i);
}
inline Outputfile::operator bool() const {
    return good();
}
inline bool Outputfile::allocate() const {
    return allocate_ofs();
}
template <typename T> const
Outputfile& Outputfile::operator<<(T& holder) const {
    test(good(),"operator<<");
    *ofs << holder;
    return *this;
}
inline ofstream& Outputfile::stream() const {
    test(good(),"stream");
    return *ofs;
}
inline bool Outputfile::allocate_ofs() const {
    return static_cast<bool>(ofs);
}
inline bool Outputfile::good() const {
    return (allocate_ofs() ? ofs->good() : false);
}

//extern Input INPUT;

#endif
