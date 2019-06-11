#ifndef CELL_TEMP_H
#define CELL_TEMP_H

#include "Cell.h"

class cell_temp1: public cell{
public:
    virtual std::istream& read_atoms(std::istream&);
};

#endif
