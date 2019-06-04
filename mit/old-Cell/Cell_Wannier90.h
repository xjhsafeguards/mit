#ifndef CELL_WANNIER90_H
#define CELL_WANNIER90_H

#include "Cell.h"

class cell_wannier90: public cell{
public:
    virtual std::istream& read(std::istream&);
};

#endif
