#include "cdftcell.h"

Inputfile& operator>>(Inputfile& inf,DFTCell& cel){
    if(inf.type()=="cif") cel.read_cif(inf.stream());
    else
        cel.test(0,"dftcell::operator>>::wrong_inputfile_type:" + inf.type());
    return inf;
}
Outputfile& operator<<(Outputfile& of, const DFTCell& cel){
    if(of.type() == "in") cel.write_in(of.stream());
    else if(of.type() == "POSCAR") cel.write_POSCAR(of.stream());
    else
        cel.test(0,"Cell::operator<<::wrong_outputfile_type:" + of.type());
    return of;
}
Outputfile& operator<<(Outputfile& of,DFTCell& cel){
    if(of.type() == "in") cel.write_in(of.stream());
    else if(of.type() == "POSCAR") cel.write_POSCAR(of.stream());
    else
        cel.test(0,"Cell::operator<<::wrong_outputfile_type:" + of.type());
    return of;
}

void DFTCell::Routine(){
    Help::help("dftcell",INPUT.type);
    switch (INPUT.type) {
        case 0:
            DFTCell::File_type_convert();
            break;
        default:
            break;
    }
}
void DFTCell::File_type_convert(){
    for( auto& file : INPUT.in_files){
        DFTCell cel;
        file >> cel;
        Outputfile of(file.name_only()+"."+INPUT.out_files_type);
        of << cel;
    }
}
