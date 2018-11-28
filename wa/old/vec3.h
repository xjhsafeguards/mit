#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>

using namespace std;

template <class T>
class Vector3
{
public:
    
    T x;
    T y;
    T z;
    
    Vector3(const T& x1=0, const T& x2=0, const T& x3=0)
    {   x = x1;y = x2;z = x3;   }
    
    Vector3(const Vector3<T> &obj)
    {   x = obj.x; y = obj.y; z = obj.z;    }
    
    void set(const T &x1,const T &x2, const T &x3)
    {   x = x1;y = x2;z = x3;   }
    
    T norm() const
    {   return sqrt(x*x+y*y+z*z);   }
    
    Vector3 uniform() const
    {   return *this/this->norm();}
    
    Vector3 operator- (const Vector3 &v2) const
    {   return Vector3(this->x-v2.x,this->y-v2.y,this->z-v2.z);   }
    
    void operator-= (const Vector3 &v2)
    {   this->set(this->x-v2.x,this->y-v2.y,this->z-v2.z);}
    
    Vector3 operator+ (const Vector3 &v2) const
    {   return Vector3(this->x+v2.x,this->y+v2.y,this->z+v2.z);   }
    
    void operator+= (const Vector3 &v2)
    {   this->set(this->x+v2.x,this->y+v2.y,this->z+v2.z);}
    
    T operator* (const Vector3 &v2) const
    // define dot
    {   return (this->x*v2.x+this->y*v2.y+this->z*v2.z); }
    
    Vector3 operator* (const double d) const
    {   return Vector3(this->x*d,this->y*d,this->z*d);}
    
    void operator*= (const double d)
    {   this->set(this->x*d,this->y*d,this->z*d);}
    
    Vector3 operator/ (const double d) const
    {   return Vector3(this->x/d,this->y/d,this->z/d);}
    
    void operator/= (const double d)
    {   this->set(this->x/d,this->y/d,this->z/d);}
    
    Vector3 cross(const Vector3 &v2) const
    {   return Vector3(y*v2.z-z*v2.y,z*v2.x-x*v2.z,x*v2.y-y*v2.x);}
    
    bool operator== (const Vector3 &v2) const
    {
        if( x == v2.x && y == v2.y && z==v2.z)
            return true;
        else
            return false;
    }
    
    T distance(const Vector3 &v2) const
    {   return (*this-v2).norm();    }
    
    T distance_BC(const Vector3 &v1,const Vector3 &BC) const //calculation distance between *this and v2 using boundary condition v3
    {
        Vector3<T> v2 = v1;
        while((this->x - v2.x) >= BC.x/2)
            v2.x += BC.x;
        while((this->x - v2.x) < -BC.x/2)
            v2.x -= BC.x;
        while((this->y - v2.y) >= BC.y/2)
            v2.y += BC.y;
        while((this->y - v2.y) < -BC.y/2)
            v2.y -= BC.y;
        while((this->z - v2.z) >= BC.z/2)
            v2.z += BC.z;
        while((this->z - v2.z) < -BC.z/2)
            v2.z -= BC.z;
        return (*this-v2).norm();
    }
    
    Vector3 shortest(const Vector3 &v1,const Vector3 &BC) const// given the shortest vector this-v1(from v1 to this)
    {
        Vector3<T> v2 = v1;
        while((this->x - v2.x) >= BC.x/2)
            v2.x += BC.x;
        while((this->x - v2.x) < -BC.x/2)
            v2.x -= BC.x;
        while((this->y - v2.y) >= BC.y/2)
            v2.y += BC.y;
        while((this->y - v2.y) < -BC.y/2)
            v2.y -= BC.y;
        while((this->z - v2.z) >= BC.z/2)
            v2.z += BC.z;
        while((this->z - v2.z) < -BC.z/2)
            v2.z -= BC.z;
        return (*this-v2);
    }
    
    void normal_BC(const Vector3 &BC)
    {
        while(x >= BC.x)
            x -= BC.x;
        while(x < 0)
            x += BC.x;
        while(y >= BC.y)
            y -= BC.y;
        while(y < 0)
            y += BC.y;
        while(z >= BC.z)
            z -= BC.z;
        while(z < 0)
            z += BC.z;
    }
    
    double angle(const Vector3 &v2,const Vector3 &v3,const Vector3 &BC) const
    // calculate the angle between v2-this-v3, return the anlge in degree
    {
        Vector3<T> v1 = v2;
        T d12 = distance_BC(v2,BC);
        T d13 = distance_BC(v3,BC);
        T d23 = v1.distance_BC(v3,BC);
        double angle = acos((d12*d12 + d13*d13 - d23*d23)/(2*d12*d13))/3.1415926535897*180;
        return angle;
    }
    
    //add 2018.3.22
    T& operator[](const int i)
    {   switch(i){
        case 0: return x; case 1: return y; case 2: return z;
        default: cerr << "Vector3 has only 3 values! " << endl; assert(0); return x;
    }}
    
    const T& operator[](const int i) const
    {   switch(i){
        case 0: return x; case 1: return y; case 2: return z;
        default: cerr << "Vector3 has only 3 values! " << endl; assert(0); return x;
    }}
    
    void operator= (const Vector3 &v2)
    {   this->set(v2.x,v2.y,v2.z);}
    
    //old functions
    
    void print(ofstream &ofs, bool endline = true) const
    {
        assert(ofs.good());
        ofs << " " << x << "  " << y << "  " << z ;
        if(endline)
            ofs << endl;
    }
    
    void print2screen(void) const
    {
        cout << "(" << x << "," << y << "," << z << ")" << endl;
    }
    
};

template <class T>
istream& operator>>(istream& is, Vector3<T> &vec)
{
    assert(is.good());
    is >> vec.x >> vec.y >> vec.z;
    return is;
}

template <class T>
ostream& operator<<(ostream& os, const Vector3<T> &vec)
{
    assert(os.good());
    //os << setw(20) << vec.x << setw(20) << vec.y << setw(20) << vec.z;
    os << vec.x << "  " << vec.y << "  " << vec.z;
    return os;
}


#endif

