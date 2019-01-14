#ifndef MD_FILEMANAGER_H
#define MD_FILEMANAGER_H


#include "Cell.h"

class MD_filemanager{
public:
    virtual void read_infile() = 0;
    virtual Cell& read_next() = 0;
    
};
#endif
