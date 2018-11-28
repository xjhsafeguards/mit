#ifndef GLOG_H
#define GLOG_H

#include "gheader.h"
#include "gmpi.h"

class Glog{
    ofstream ofs_log;
    clock_t start_time,end_time;
public:
    Glog(){}

    void init();
    void end();
    
    //set
    ofstream& setprecision(size_t p=10);
    
    ostream& out_setprecision(size_t p=10);
    
    //write
    Glog* write(const string& s);
    Glog* write_line(const string& s);
    ofstream& stream();
    
    ostream& out_stream();
};

extern Glog *GLOG;

#endif
