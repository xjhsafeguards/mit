#ifndef CELL_H
#define CELL_H

#include <string>
#include <vector>
#include <map>

#include "Math_linearalgebra.h"


class box;
class position;
class frac_position;
class cart_position;

class Cell;

class box{
public:
    typedef Matrix3<double>     data_type;
    typedef Vector3<double>     diagonal_type;
    typedef Vector3<double>     pos_type;
    
    //constructor
    box(const data_type& in_cell): celldm(in_cell), inverse_celldm(in_cell.inverse()), is_orthorhombic(in_cell.orthorhombic()){}
    box(const box&) = default;
    box(box&&) = default;
    box& operator=(const box&) = default;
    box& operator=(box&&) = default;
    
    //operatrions
    const data_type& data() const{
        return celldm;
    }
    const data_type& inverse() const{
        return inverse_celldm;
    }
    diagonal_type diagonal() const{
        assert(is_orthorhombic);
        return celldm.diagonal();
    }
    bool orthorhombic() const{
        return is_orthorhombic;
    }
    
protected:
    data_type celldm;
    data_type inverse_celldm;
    bool      is_orthorhombic;
};

class position{
public:
    typedef Vector3<double>         data_type;
    typedef std::shared_ptr<box>    box_ptr;
    
    //constructor
    //position(){};
    position(const data_type& in_p,const box_ptr& in_b): pos(in_p),box(in_b){}
    position(const position&) = default;
    position(position&&) = default;
    position& operator=(const position&) = default;
    position& operator=(position&&) = default;
    
    //operations

    //output position in fraction or cartesian
    virtual data_type frac() const = 0;
    virtual data_type cart() const = 0;
    // distance between this and p
    virtual double distance(const position& p) const = 0;
    // angle between p1-this-p2 in radian
    virtual double angle(const position& p1, const position& p2) const = 0;
    
    
protected:
    data_type pos;
    box_ptr box;
    
    const data_type& data() const{
        return pos;
    }
    // compare if this and p in the the same box
    bool same_box(const position& p) const{
        return p.box == box;
    }
    data_type cartesian(const data_type& p) const{
        return p*box->data();
    }
    data_type fraction(const data_type& p) const{
        return p*box->inverse();
    }
    
    //folder operatrions
};



class frac_position : public position{
public:
    frac_position(const data_type& in_p,const box_ptr& in_b,bool is_f=true): position(in_p,in_b){
        if(!is_f)
            pos = fraction(pos);
            }
    //frac_position(const frac_position&) = default;
    //frac_position(frac_position&&) = default;
    frac_position& operator=(const frac_position&) = default;
    frac_position& operator=(frac_position&&) = default;
    
    virtual data_type frac() const{
        return pos;
    }
    virtual data_type cart() const{
        return cartesian(pos);
    }
    
    // distance between this and p
    virtual double distance(const position& p) const;
    // angle between p1-this-p2 in degree
    virtual double angle(const position& p1, const position& p2) const;

protected:
};

class cart_position: public position{
public:
    using position::position;
    cart_position(const data_type& in_p,const box_ptr& in_b,bool is_c=true): position(in_p,in_b){
        if(!is_c)
            pos = cartesian(pos);
    }
    cart_position& operator=(const cart_position&) = default;
    cart_position& operator=(cart_position&&) = default;
    
    virtual data_type frac() const{
        return fraction(pos);
    }
    virtual data_type cart() const{
        return pos;
    }
    
    // distance between this and p
    virtual double distance(const position& p) const;
    // angle between p1-this-p2 in degree
    virtual double angle(const position& p1, const position& p2) const;
    
protected:
};


class Cell{
public:
    typedef std::vector<std::shared_ptr<position> > pos_folder;
    
    //get a pos_folder with certain character
    const pos_folder& get_folder(std::string character){
        return data[character];
    }
    
protected:
    std::shared_ptr<box> b;
    std::map<std::string,std::vector<std::shared_ptr<position> > > data;
    
    pos_folder& operator[](std::string character){
        return data[character];
    }
};

#endif

