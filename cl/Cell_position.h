#ifndef CELL_POSITION_H
#define CELL_POSITION_H

#include <string>
#include <vector>
#include <memory>
#include <iostream>


#include "Math_linearalgebra.h"


class box;
class position;
class frac_position;
class cart_position;

//class Cell;

class box{
public:
    typedef Matrix3<double>     data_type;
    typedef Vector3<double>     diagonal_type;
    typedef Vector3<double>     pos_type;
    
    //constructor
    box(const data_type& in_cell): celldm(in_cell), inverse_celldm(in_cell.inverse()), is_orthorhombic(in_cell.orthorhombic()){}
    box(std::istream& is){
        is >> celldm;
        inverse_celldm = celldm.inverse();
        is_orthorhombic = celldm.orthorhombic();
    }
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
    
    //set
    void set_unit(double in_unit){
        celldm = celldm*in_unit;
    }
    
    //read
    double volume() const{
        return celldm[0].cross(celldm[1])*celldm[2];
    }
    
    std::ostream& write(std::ostream& os){
        os << celldm;
        return os;
    }
    std::ostream& write_inverse(std::ostream& os){
        os << inverse_celldm;
        return os;
    }
    
protected:
    data_type celldm;
    data_type inverse_celldm;
    bool      is_orthorhombic;
};

class position{
public:
    typedef std::string             type_type;
    typedef int                     label_type;
    typedef Vector3<double>         data_type;
    typedef std::shared_ptr<box>    box_ptr;

    //constructor
    position(){};
    position(const data_type& in_p,const box_ptr& in_b,const std::string& in_type="", int in_label=-1): pos(in_p),boxp(in_b),type(in_type),label(in_label){}
    position(const position&) = default;
    position(position&&) = default;
    position& operator=(const position&) = default;
    position& operator=(position&&) = default;
    
    //operations

    //set operations
    inline void set_box_ptr(const box_ptr& in_b);
    //inline void set_box_ptr(std::istream&);
    inline void set_position(double x,double y,double z);
    inline void set_position(std::istream&);
    inline void set_type(const type_type& in_type);
    inline void set_type(std::istream&);
    inline void set_label(const label_type& in_label);
    inline void set_unit(double in_unit);
    
    //check_type
    inline bool check_type(const std::string& in_type) const;
    
    //get
    inline const type_type& get_type() const;
    
    //output position in fraction or cartesian
    virtual data_type frac() const = 0;
    virtual data_type cart() const = 0;
    // distance between this and p
    virtual double distance(const position& p) const = 0;
    // angle between p1-this-p2 in radian
    virtual double angle(const position& p1, const position& p2) const = 0;
    
    // stream operations
    virtual std::ostream& write(std::ostream& os){
        os << pos;
        return os;
    }
    virtual std::ostream& write_frac(std::ostream& os){
        os << frac();
        return os;
    }
    virtual std::ostream& write_cart(std::ostream& os){
        os << cart();
        return os;
    }
    
protected:
    type_type type;
    label_type label;
    data_type pos;
    box_ptr boxp;
    
    const data_type& data() const{
        return pos;
    }
    // compare if this and p in the the same box
    bool same_box(const position& p) const{
        return p.boxp == boxp;
    }
    data_type cartesian(const data_type& p) const{
        return p*boxp->data();
    }
    data_type fraction(const data_type& p) const{
        return p*boxp->inverse();
    }
    
};



class frac_position : public position{
public:
    using position::position;
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


//set operations
inline void position::set_box_ptr(const box_ptr& in_b){
    boxp = in_b;
}
inline void position::set_position(double x,double y,double z){
    pos = data_type(x,y,z);
}
inline void position::set_position(std::istream& is){
    is >> pos;
}
inline void position::set_type(const type_type& in_type){
    type = in_type;
}
inline void position::set_type(std::istream& is){
    is >> type;
}
inline void position::set_label(const label_type& in_label){
    label = in_label;
}
inline void position::set_unit(double in_unit){
    pos *= in_unit;
}

//check_type
inline bool position::check_type(const std::string& in_type) const{
    return in_type == type;
}

inline const std::string& position::get_type() const{
    return type;
}


/*
class Cell{
public:
    typedef std::vector<std::shared_ptr<position> > pos_folder
    
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
 */

#endif

