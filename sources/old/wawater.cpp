#include "wawater.h"
//Waters
Outputfile& operator<<(Outputfile& of,Waters& wa){
    if(of.type() == "txt") wa.write_waters(of.stream());
    else
        wa.test(0,"operator<<::wrong_outputfile_type:" + of.type());
    return of;
}
Waters::Waters(Cell& in_cel,bool ifsave) {
    ifsave ? cel = in_cel.save() : cel = make_shared<Cell>(in_cel);
}
void Waters::set_waters(){
    int O,H;
    try{
        O = cel->aindex("O");
        H = cel->aindex("H");
    }catch(runtime_error){
        test(0,"set_waters::need_O_and_H");
    }
    waters.clear();
    int nO = cel->anum("O"), nH = cel->anum("H");
    for(int i=0; i<nO; i++){
        Water tmpwater(cel,i);
        tmpwater.set_water();
        waters.push_back(make_shared<Water>(tmpwater));
        if(tmpwater.nH()!=2){
            ions.push_back(i);
        }
    }
}
void Waters::write_waters(ostream& os){
    os << cel->snapshot() << " " << cel->time() << endl;
    for( const auto& wa: waters){
        wa->write_index(os);
        //os << "   " << wa->posO();
        //for( int i=0; i<wa->nH(); i++)
        //    os << " " << wa->posH(i);
        os << endl;
    }
}
void Waters::test(bool condition,string information) const{
    if(!condition)
        throw(runtime_error("Waters::Fail_"+information));
}
