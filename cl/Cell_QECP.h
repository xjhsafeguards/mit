#ifndef CELL_QECP_H
#define CELL_QECP_H

#include "Cell.h"

class cell_qecp: public cell{
public:
    virtual std::istream& read_atoms(std::istream&);
};

#endif

