#include "gfile.h"

File::File(string name) : filename(name){
    set_suffix_filetype();
}
File::File(string name,string type) : filename(name),filetype(type){
}
string File::name_only() const{
#ifndef N_TEST
    test(allocate_name(),"name_only");
#endif
    auto dot = find(filename.crbegin(),filename.crend(),'.');
    if(dot!=filename.crend())
        return string(filename.cbegin(),dot.base()-1);
    else
        return name();
}
void File::set_suffix_filetype()
{
#ifndef N_TEST
    test(allocate_name(),"set_suffix_filetype");
#endif
    auto dot = find(filename.crbegin(),filename.crend(),'.');
    filetype = string(dot.base(),filename.cend());
}
istream& operator>>(istream& is, File& file){
    is >> file.filename;
    file.open();
    return is;
}
ostream& operator<<(ostream& os, File& file){
    os << file.filename;
    return os;
}

Inputfile::Inputfile(string name): File(name),ifs(make_shared<ifstream>(name)) {
    test(good(),"Open " + filename + " to read");
    log("open file " + filename + " to read");
}
Inputfile::Inputfile(string name,string type): File(name,type),ifs(make_shared<ifstream>(name)) {
    test(good(),"Open " + filename + " to read");
    log("open file " + filename + " to read");
}
Inputfile::~Inputfile() {
    if(allocate_ifs() and ifs.unique())
        log("close file " + filename);
}
void Inputfile::open()
{
    ifs = make_shared<ifstream>(filename);
    test(good(),"Open " + filename + " to read");
    log("open file " + filename + " to read");
    set_suffix_filetype();
}
void Inputfile::open(string name)
{
    filename = name;
    open();
}
void Inputfile::seek_top(bool show_log){
    ifs->seekg(0,ifs->beg);
    if(show_log)
        log(" file " + filename + " seek_top");
}
Outputfile::Outputfile(string name): File(name),ofs(make_shared<ofstream>(name)) {
    test(good(),"Open " + filename + " to write");
    log("open file " + filename + " to write");
    setprecision();
}
Outputfile::Outputfile(string name,string type): File(name,type),ofs(make_shared<ofstream>(name)) {
    test(good(),"Open " + filename + " to write");
    log("open file " + filename + " to write");
    setprecision();
}
Outputfile::~Outputfile() {
    if(allocate_ofs() and ofs.unique())
        log("close file " + filename);
}
void Outputfile::open()
{
    ofs = make_shared<ofstream>(filename);
    test(good(),"Open " + filename + " to write");
    log("open file " + filename + " to write");
    setprecision();
    set_suffix_filetype();
}
void Outputfile::open(string name)
{
    filename = name;
    open();
}

