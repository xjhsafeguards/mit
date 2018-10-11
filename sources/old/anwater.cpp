#include "anwater.h"
#include "input.h"

double Water::OH_distance = 1.26;

Water::Water(shared_ptr<Cell> in_cel): cel(in_cel){
}
Water::Water(shared_ptr<Cell> in_cel,  int in_indexO): cel(in_cel),indexO(in_indexO){
}
void Water::write_index(ostream& os) const{
    os << setw(4) << indexO;
    for( auto& H : indexH)
        os << setw(4) << H;
}
void Water::set_water(){
    int O = cel->aindex("O");
    int H = cel->aindex("H");
    for(int j=0; j< cel->anum("H"); j++){
        if(cel->distance(O,indexO,H,j) < OH_distance){
            indexH.push_back(j);
        }
    }
}
void Water::test(bool condition,string information) const{
    if(!condition)
        throw(runtime_error("Water::Fail_"+information));
}

