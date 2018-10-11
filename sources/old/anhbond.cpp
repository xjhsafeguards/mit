#include "anhbond.h"

bool Hbond::check_accept(int i) const{
    auto result = find(index_accept.begin(),index_accept.end(),i);
    return result!=index_accept.end();
}
bool Hbond::check_donate(int i) const{
    auto result = find(index_donate.begin(),index_donate.end(),i);
    return result!=index_donate.end();
}
bool Hbond::check_bond(int i) const{
    return check_accept(i) or check_donate(i);
}

void Hbond::write_detail_header(ostream& os){
    os << setw(5) << "nA" << setw(5) << "nD" << setw(5) << "A1" << setw(5) << "A2" << setw(5) << "A3" << setw(5) << "A4" << setw(5) << "A5" << setw(5) << "A6" << setw(5) << "A7" << setw(5) << "A8" << setw(5) << "D1" << setw(5) << "D2" << setw(5) << "D3" << setw(5) << "D4" << endl;
}
void Hbond::write_detail(ostream& os) const{
    os << setw(5) << naccept() << setw(5) << ndonate();
    for(int j=0; j<8; j++){
        if(j<naccept())
            os << setw(5) << index_accept[j];
        else
            os << setw(5) << "-1";
    }
    for(int j=0; j<4; j++){
        if(j<ndonate())
            os << setw(5) << index_donate[j];
        else
            os << setw(5) << "-1";
    }
    os << endl;
}

void Hbond::test(bool condition,string information) const{
    if(!condition)
        throw(runtime_error("Hbond::Fail_"+information));
}
