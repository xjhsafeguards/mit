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

#include "vec3.h"
#include "mat3.h"
#include "math.h"
#include "phy.h"


using namespace std;

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

inline void Skip_lines(ifstream& ifs, int n)
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

inline void my_assert(bool test, string output="Some thing went wrong!")
{
    if(!test)
    {
        cerr << output << endl;
        exit(1);
    }
}

inline void file_assert(ifstream& ifs,string file_name,int type = 0)
{
    if(ifs.fail())
    {
        if(type == 0)
            cerr << "Can not file: " << file_name  << " to read!" << endl;
        else
            cerr << file_name << endl;
        exit(1);
    }
}

inline void file_assert(ofstream& ofs,string file_name,int type = 0)
{
    if(ofs.fail())
    {
        if(type == 0)
            cerr << "Can not file: " << file_name  << " to write!" << endl;
        else
            cerr << file_name << endl;
        exit(1);
    }
}









#endif
