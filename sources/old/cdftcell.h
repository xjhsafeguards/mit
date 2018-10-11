#ifndef CDFTCELL_H
#define CDFTCELL_H

#include "ccell.h"

class DFTCell;
Inputfile& operator>>(Inputfile&, DFTCell&);
Outputfile& operator<<(Outputfile&, const DFTCell&);
Outputfile& operator<<(Outputfile&,DFTCell&);

class DFTCell: public Cell{
    friend Inputfile& operator>>(Inputfile&, DFTCell&);
    friend Outputfile& operator<<(Outputfile&, const DFTCell&);
    friend Outputfile& operator<<(Outputfile&,DFTCell&);
private:
    
public:
    static void Routine();
    static void File_type_convert();
    
    DFTCell() = default;
    
    
    
    //2018.6.14 add by shuo
    DFTCell& operator<<(Inputfile& inf){
        if(inf.type()=="cif")
            this->read_cif(inf.stream());
        return *this;
    }
    DFTCell& operator>>(ostream& os){
        this->write_in(os);
        return *this;
    }
};

#endif
