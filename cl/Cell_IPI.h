#ifndef CELL_IPI_H
#define CELL_IPI_H

#include <string>

#include "Cell.h"
#include "Physics_unit.h"

class cell_ipi: public cell{
    int cell_type;
public:
    cell_ipi(int in_type=0): cell_type(in_type){}

    virtual std::istream& read(std::istream&);
};

#endif
