#ifndef CELL_QECP_H
#define CELL_QECP_H

#include <cassert>
#include <vector>
#include <string>
#include <initializer_list>

#include "Cell.h"
#include "Physics_unit.h"

class cell_qecp: public cell{
    std::vector<int> na;
    std::vector<std::string> types;
public:
    cell_qecp(){}
    cell_qecp(std::initializer_list<int> li_na,std::initializer_list<std::string> li_types): na(li_na), types(li_types){
        assert(na.size()==types.size());
    }
    
    virtual std::istream& read_box(std::istream&);
    virtual std::istream& read_atoms(std::istream&);
};

#endif

