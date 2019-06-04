#include "Cell_position.h"

// distance between this and p
double frac_position::distance(const position& p) const{
    assert(same_box(p));
    return cartesian(pos.shortest_BC(p.frac())).norm();
}
// angle between p1-this-p2 in radian
double frac_position::angle(const position& p1, const position& p2) const{
    assert(same_box(p1));
    assert(same_box(p2));
    data_type v1 = pos.shortest_BC(p1.frac());
    data_type v2 = pos.shortest_BC(p2.frac());
    return cartesian(v1).vangle(cartesian(v2));
}

// distance between this and p
double cart_position::distance(const position& p) const{
    assert(same_box(p));
    return pos.distance_BC(p.cart(),boxp->diagonal());
}
// angle between p1-this-p2 in radian
double cart_position::angle(const position& p1, const position& p2) const{
    assert(same_box(p1));
    assert(same_box(p2));
    return pos.angle_BC(p1.cart(),p2.cart(),boxp->diagonal());
}

