#include "Cell.h"

Atom::pos_type Atom::position() const{
    return cel.atom_position(index);
}

double Atom::distance(int j) const{
    return cel.distance(index,j);
}

double Atom::angle(int j, int k) const{
    return cel.angle(index,j,k);
}

/***************** Cell *******************/

Atom Cell::atom(int index) const {return Atom(*this,index);}
Box  Cell::box() const {return Box(*this);}

double Cell::distance(int i,int j) const{
    if(is_frac)
        return cal_distance_f(atom_position(i),atom_position(j),cell_parameters);
    return cal_distance_c(atom_position(i),atom_position(j),cell_parameters);
}

double Cell::angle(int i,int j, int k) const{
    if(is_frac)
        return cal_angle_f(atom_position(i),atom_position(j),atom_position(k),cell_parameters);
    return cal_angle_c(atom_position(i),atom_position(j),atom_position(k),cell_parameters);
}

double Cell::volume() const{
    return cell_parameters.row(0).cross(cell_parameters.row(1)).dot(cell_parameters.row(2));
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
        ++iter;
    }
    return std::move(v);
}

typename Cell::pos_type Cell::shortest_vector(const pos_type& cp1,const pos_type& cp2,const box_type& bc) const{
    pos_type v = cp1 - cp2;
    auto iter=v.data();
    for(int i=0;i<3;i++){
        double bcv=bc(i,i);
        double hbcv=bcv/2;
        while(*iter<=-hbcv) *iter+=bcv;
        while(*iter>hbcv) *iter-=bcv;
        ++iter;
    }
    return std::move(v);
}
