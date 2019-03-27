#include "Cell.h"

double Atom::distance() const{
    return cel.distance();
}

/***************** Cell *******************/
double Cell::distance() const{
    return 1;
    
}

double Cell::angle() const{
    return 1;
    
}

void Cell::to_frac(){
    if(is_frac)
        return;
    else{
        positions *= cell_parameters.inverse();
        is_frac=true;
    }
    return;
}
void Cell::to_cart(){
    if(!is_frac)
        return;
    else{
        positions *= cell_parameters;
        is_frac=false;
    }
    return;
}
typename Cell::pos_type Cell::shortest_fvector(const pos_type& fp1,const pos_type& fp2) const{
    pos_type v = fp1 - fp2;
    /*
    while(v(0)<=0) v(0)+=1;
    while(v(0)>1) v(0)-=1;
    while(v(1)<=0) v(1)+=1;
    while(v(1)>1) v(1)-=1;
    while(v(2)<=0) v(2)+=1;
    while(v(2)>1) v(2)-=1;
     */
    auto iter=v.data();
    for(int i=0;i<3;i++){
        while(*iter<=-0.5) *iter+=1;
        while(*iter>0.5) *iter-=1;
    }
    return std::move(v);
}
