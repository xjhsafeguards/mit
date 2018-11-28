#ifndef GHELP_H
#define GHELP_H

#include "gfun.h"

class Help
{
public:
    
    static void Routine(string calculation);
    static void help(string calculation,int method);
    static map<string,vector<string> > dic;
    
    static void title(string calculation);
    static void content(string calculation,int method);
    static int  size(string calculation);
};

#endif
