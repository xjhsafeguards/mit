#ifndef GFUN_H
#define GFUN_H

#include <stdlib.h>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cassert>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <map>
#include <memory>
#include <typeinfo>
#include <tuple>

using namespace std;
// for mpi
extern int NCOR;
extern int RANK;
// for log
extern ofstream ofs_log;

#include "gmath.h"


template<class T>
void convstring(const string &s,T &temp)
{
    stringstream ss;
    ss << s;
    ss >> temp;
}

template<class T>
void convstring(const T &temp,string &s)
{
    stringstream ss;
    ss << temp;
    ss >> s;
}

template <class T>
void Read_value(ifstream& ifs, T& data)
{
    assert(ifs.good());
    ifs >> data;
    ifs.ignore(200,'\n');
}

inline void Seek_value(istream& ifs,string title)
{
    assert(ifs.good());
    string buffer;
    while(ifs){
        ifs >> buffer;
        if(buffer == title)
            return;
        ifs.ignore(200,'\n');
    }
    throw(runtime_error("Seek_value::no_"+title+"_found"));
}

template <class T>
void Seek_value(istream& ifs,string title,T& data)
{
    assert(ifs.good());
    string buffer;
    while(ifs){
        ifs >> buffer;
        if(buffer == title){
            ifs >> data;
            ifs.ignore(200,'\n');
            return;
        }
        ifs.ignore(200,'\n');
    }
    throw(runtime_error("Seek_value::no_"+title+"_found"));
}

template <class T>
void Read_values(ifstream& ifs, T *data, int n)
{
    assert(ifs.good());
    for(int i=0;i<n;i++)
        ifs >> data[i];
    ifs.ignore(200,'\n');
}

template <class T>
void Read_line(istream &is, T *data, int n)
{
    assert(is.good());
    string s, tmp;
    int count=0;
    getline(is,s);
    s.push_back('\t');
    for(auto ch : s)
    {
        if(isspace(ch) and !tmp.empty() )
        {
           if (count < n)
           {
               convstring(tmp,data[count++]);
               tmp.clear();
           }
           else return;
        }
        else if(!isspace(ch))
            tmp.push_back(ch);
    }
}

inline void Skip_lines(istream& ifs, int n)
//skip n lines in ifs
{
    assert(ifs.good());
    for(int i=0;i<n;i++)
        ifs.ignore(500,'\n');
}

inline void my_throw(bool test,const string what="Some thing went wrong!")
{
    if(!test)
        throw(std::runtime_error(what));
}

inline bool isspace(string str){
    bool space = true;
    for(auto it=str.cbegin(); it != str.cend(); it++)
        if(!isspace(*it))  space = false;
    return space;
}


template <typename T>
shared_ptr<T> value_copy(const shared_ptr<T> &obj){
    //T tmpobj(*obj);
    return make_shared<T>(*obj);
}

#endif
