#ifndef MATH_LINEARALGEBRA_H
#define MATH_LINEARALGEBRA_H

#include <cassert>
#include <exception>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "Math_const.h"

template<typename> class Vector3;
template <typename> class Matrix3;
template <typename T> std::istream& operator>>(std::istream& is, Vector3<T> &vec);
template <typename T> std::ostream& operator<<(std::ostream& os, const Vector3<T> &vec);
template <typename T> void swap(Vector3<T>& v1,Vector3<T>& v2);
template <typename T> std::ostream& operator<<(std::ostream&, const Matrix3<T>&);
template <typename T> std::istream& operator>>(std::istream&, Matrix3<T>&);
template <typename T> typename Matrix3<T>::row_type operator*(const Vector3<T>&, const Matrix3<T>&);
template <typename T> void swap(Matrix3<T>& m1,Matrix3<T>& m2);

/* ***********************************************************
 epsilon
 *************************************************************/
inline int epsilon(const int i, const int j, const int k)
{
    assert(i>0 and i<4);
    assert(j>0 and j<4);
    assert(k>0 and k<4);
    if(i==j or i==k or k==j)
        return 0;
    else
        return -(i-j)*(j-k)*(i-k)/2;
}

/* ***********************************************************
    Vector3
*************************************************************/
template <typename T>
class Vector3
{
    friend std::istream& operator>> <T> (std::istream&, Vector3&);
    friend std::ostream& operator<< <T> (std::ostream&, const Vector3&);
    friend class Matrix3<T>;
    friend void swap<T>(Vector3& v1,Vector3 &v2);
    
    T x;
    T y;
    T z;
    
public:
    //type
    typedef T value_type;
    //constructor
    Vector3(const T& x1=0, const T& x2=0, const T& x3=0): x(x1),y(x2),z(x3){}
    Vector3(const Vector3<T> &obj): x(obj.x),y(obj.y),z(obj.z){}
    Vector3(Vector3<T> &&obj) noexcept: x(std::move(obj.x)),y(std::move(obj.y)),z(std::move(obj.z)){}
    ~Vector3() = default;
    
    Vector3& operator= (Vector3 v2) &
    {
        swap(*this,v2);
        return *this;
    }
    //assignment
    void set(const T &x1,const T &x2, const T &x3)
    {
        x = x1;
        y = x2;
        z = x3;
    }
    //read
    T& operator[](const int i)
    {
        switch(i){
            case 0: return x;
            case 1: return y;
            case 2: return z;
        }
        vthrow("[]::out_of_range");
        return x;
    }
    const T& operator[](const int i) const
    {
        switch(i){
            case 0: return x;
            case 1: return y;
            case 2: return z;
        }
        vthrow("[]::out_of_range");
        return x;
    }
    //algorithm
    Vector3 operator+ (const Vector3 &v2) const
    {
        return Vector3(this->x+v2.x,this->y+v2.y,this->z+v2.z);
    }
    Vector3 operator- (const Vector3 &v2) const
    {
        return Vector3(this->x-v2.x,this->y-v2.y,this->z-v2.z);
    }
    Vector3 operator* (const double d) const
    {
        return Vector3(this->x*d,this->y*d,this->z*d);
    }
    Vector3 operator/ (const double d) const
    {
        return Vector3(this->x/d,this->y/d,this->z/d);
    }
    Vector3& operator+= (const Vector3 &v2)
    {
        this->set(this->x+v2.x,this->y+v2.y,this->z+v2.z);
        return *this;
    }
    Vector3& operator-= (const Vector3 &v2)
    {
        this->set(this->x-v2.x,this->y-v2.y,this->z-v2.z);
        return *this;
    }
    Vector3& operator*= (const double d)
    {
        this->set(this->x*d,this->y*d,this->z*d);
        return *this;
    }
    Vector3& operator/= (const double d)
    {
        this->set(this->x/d,this->y/d,this->z/d);
        return *this;
    }
    //compare
    bool operator== (const Vector3 &v2) const
    {
        if( x == v2.x && y == v2.y && z==v2.z)
            return true;
        else
            return false;
    }
    //other operations
    T norm() const
    {
        return sqrt(x*x+y*y+z*z);
    }
    Vector3 uniform()
    {   return *this/this->norm();}
    //dot
    T operator* (const Vector3 &v2) const
    {
        return (this->x*v2.x+this->y*v2.y+this->z*v2.z);
    }
    T dot(const Vector3 &v2) const
    {
        return (this->x*v2.x+this->y*v2.y+this->z*v2.z);
    }
    Vector3 cross(const Vector3 &v2) const
    {
        return Vector3(y*v2.z-z*v2.y,z*v2.x-x*v2.z,x*v2.y-y*v2.x);
    }
    // calculate the angle between this-v2 , return the anlge in degree
    double vangle(const Vector3 &v2) const
    {
        return acos(this->dot(v2)/this->norm()/v2.norm())/PI*180;
    }
    
    //point like actions
    T distance(const Vector3 &v2) const
    {
        return (*this-v2).norm();
    }
    
//calculate with cuboid boxs whose diagonal element is BC
    //put this under boudary condition
    void normal_BC(const Vector3 &BC=Vector3(1,1,1))
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
    //point like actions
    //calculation distance between *this and v1 using boundary condition BC
    T distance_BC(const Vector3 &v1,const Vector3 &BC=Vector3(1,1,1)) const
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
    // calculate the shortest vector this-v1(from v1 to this)
    Vector3 shortest_BC(const Vector3 &v1,const Vector3 &BC=Vector3(1,1,1)) const
    {
        Vector3 v2 = v1;
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
    double angle_BC(const Vector3 &v2,const Vector3 &v3,const Vector3 &BC=Vector3(1,1,1)) const
    // calculate the angle between v2-this-v3, return the anlge in degree
    {
        Vector3 v21 = shortest_BC(v2,BC);
        Vector3 v31 = shortest_BC(v3,BC);
        return v31.vangle(v21);
    }
    double angle_BC2(const Vector3 &v2,const Vector3 &v3,const Vector3 &BC=Vector3(1,1,1)) const
    // calculate the angle between v2-this-v3, return the anlge in degree
    {
        Vector3 v1 = v2;
        T d12 = distance_BC(v2,BC);
        T d13 = distance_BC(v3,BC);
        T d23 = v1.distance_BC(v3,BC);
        return acos((d12*d12 + d13*d13 - d23*d23)/(2*d12*d13))/PI*180;
    }
    void print(std::ofstream &ofs, bool endline = true) const
    {
        assert(ofs.good());
        ofs << " " << x << "  " << y << "  " << z ;
        if(endline)
            ofs << std::endl;
    }
    
    void print2screen(void) const
    {
        std::cout << "(" << x << "," << y << "," << z << ")" << std::endl;
    }
    
private:
    void vthrow(std::string information) const
    {
        throw(std::runtime_error("VECTOR3::" + information));
    }
    
};

template <typename T>
std::istream& operator>>(std::istream& is, Vector3<T> &vec)
{
    is >> vec.x >> vec.y >> vec.z;
    return is;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vector3<T> &vec)
{
    os << std::setw(os.precision()+4) << vec.x << " " << std::setw(os.precision()+4) << vec.y << " " << std::setw(os.precision()+4) << vec.z;
    return os;
}
template <typename T>
void swap(Vector3<T>& v1,Vector3<T>& v2)
{
    using std::swap;
    swap(v1.x,v2.x);
    swap(v1.y,v2.y);
    swap(v1.z,v2.z);
}
/* ***********************************************************
 End of Vector3
*************************************************************/

/* ***********************************************************
 Matrix3
 *************************************************************/
template <class T>
class Matrix3
{
public:
    //type
    typedef T value_type;
    typedef Vector3<T> row_type;
    typedef Vector3<T> col_type;
    
    friend std::ostream& operator<< <T> (std::ostream&, const Matrix3&);
    friend std::istream& operator>> <T> (std::istream&, Matrix3&);
    friend typename Matrix3<T>::row_type operator* <T>(const Vector3<T>&, const Matrix3&);
    friend void swap<T>(Matrix3 &m1,Matrix3& m2);
    
private:
    Vector3<Vector3<T> > Mat;
    
public:

    //constructor
    Matrix3(const T &x11=0,const T &x12=0,const T &x13=0,
            const T &x21=0,const T &x22=0,const T &x23=0,
            const T &x31=0,const T &x32=0,const T &x33=0): Mat(Vector3<T>(x11,x12,x13),Vector3<T>(x21,x22,x23),Vector3<T>(x31,x32,x33)) {}
    Matrix3(const Vector3<T>& l1, const Vector3<T>& l2, const Vector3<T>& l3): Mat(l1,l2,l3) {}
    Matrix3(const Vector3<Vector3<T> > &mat): Mat(mat) {}
    Matrix3(const Matrix3<T> &obj): Mat(obj.Mat) {}
    Matrix3(Matrix3&& obj) noexcept: Mat(std::move(obj.Mat)){}
    
    Matrix3& operator= (Matrix3<T> m2) &
    {
        swap(*this,m2);
        return *this;
    }
    //assignment
    Matrix3& set(const T &x11=0,const T &x12=0,const T &x13=0,
                const T &x21=0,const T &x22=0,const T &x23=0,
                const T &x31=0,const T &x32=0,const T &x33=0)
    {
        Mat = Vector3<Vector3<T> >(Vector3<T>(x11,x12,x13),Vector3<T>(x21,x22,x23),Vector3<T>(x31,x32,x33));
        return *this;
    }
    //Matrix3& operator= (const Matrix3<T>& m2)
    //{
    //    Mat = m2.Mat;
    //    return *this;
    //}
    //Read
    row_type& operator[] (const int& i) //return i+1th row i(0,1,2)
    {
        return Mat[i];
    }
    const row_type& operator[] (const int& i) const    //return i+1th row i(0,1,2)
    {
        return Mat[i];
    }
    row_type row (const int& i, const bool mod= false) const    //return ith row i(1,2,3)
    {
        return Mat[i-1];
    }
    col_type col (const int& i) const    //return ith column i(1,2,3)
    {
        return Vector3<T>(Mat[0][i-1],Mat[1][i-1],Mat[2][i-1]);
    }
    value_type value(const int &i, const int &j,const bool mod= false) const    //return Mat(i,j) i,j(1,2,3)
    {
        return Mat[i-1][j-1];
    }
    //arithmetic
    Matrix3 operator+ (const Matrix3<T> &m2)
    {
        return Mat + m2.Mat;
    }
    Matrix3 operator- (const Matrix3<T> &m2)
    {
        return Mat - m2.Mat;
    }
    Matrix3 operator* (const double& num) const
    {
        Matrix3<T> tmp=*this;
        for(int i=0;i<3;++i)
            tmp[i]*=num;
        return std::move(tmp);
    }
    Matrix3 operator/ (const double& num) const
    {
        Matrix3<T> tmp=*this;
        for(int i=0;i<3;++i)
            tmp[i]/=num;
        return std::move(tmp);
    }
    col_type operator* (const row_type& v) const
    {
        return col_type(Mat[0]*v,Mat[1]*v,Mat[2]*v);
    }
    Matrix3 operator* (const Matrix3<T>& m2) const
    {   return Matrix3<T>(Vector3<T>(Mat[0]*m2.col(1),Mat[0]*m2.col(2),Mat[0]*m2.col(3)),
                          Vector3<T>(Mat[1]*m2.col(1),Mat[1]*m2.col(2),Mat[1]*m2.col(3)),
                          Vector3<T>(Mat[2]*m2.col(1),Mat[2]*m2.col(2),Mat[2]*m2.col(3)));
    }
    T det() const
    {
        T result=0;
        for(int i=1;i<=3;++i)
            for(int j=1;j<=3;++j)
                for(int k=1;k<=3;k++)
                    result += epsilon(i,j,k)*value(1,i)*value(2,j)*value(3,k);
        return std::move(result);
    }
    T cofact(const int &a, const int &b) const;    //return Cij i,j(1,2,3)
    
    //operations
    Matrix3<T> transposition() const
    {
        return Matrix3<T>(col(1),col(2),col(3));
    }
    Matrix3<T> inverse() const
    {
        assert( det()!=0 );
        Matrix3<T> result(cofact(1,1),cofact(1,2),cofact(1,3),
                          cofact(2,1),cofact(2,2),cofact(2,3),
                          cofact(3,1),cofact(3,2),cofact(3,3));
        return result.transposition()/det();
    }
    Vector3<T> diagonal() const{
        return Vector3<T>(Mat[1][1],Mat[2][2],Mat[3][3]);
    }
    bool orthorhombic() const{
        for(int i=1;i<=3;i++)
            for(int j=1;j<=3;j++)
                if(i!=j)
                    if( value(i,j) != 0)
                        return false;
        return true;
    }
    
    //old
    T cofact2(const int &i, const int &j) const     //return Cij i,j(1,2,3)
    {
        std::vector<T> d4;
        for(int x=1; x<=3; x++)
            for(int y=1; y<=3; y++)
                if( x != i and y != j)
                    d4.push_back(value(x,y));
        assert(d4.size()==4);
        return std::pow(-1,i+j)*(d4[0]*d4[3]-d4[1]*d4[2]);
    }
private:
    void mthrow(std::string information) const
    {
        throw(std::runtime_error("Matrix3::" + information));
    }
};

template <class T>
std::istream& operator>>(std::istream &is,Matrix3<T> &m)
{
    is >> m[0] >> m[1] >> m[2];
    return is;
}

template <typename T>
T Matrix3<T>::cofact(const int &a, const int &b) const     //return Cij i,j(1,2,3)
{
    T result=0;
    switch (a) {
        case 1:
            switch (b) {
                case 1:
                    result = value(2,2)*value(3,3) - value(2,3)*value(3,2);
                    break;
                case 2:
                    result = value(2,3)*value(3,1) - value(2,1)*value(3,3);
                    break;
                case 3:
                    result = value(2,1)*value(3,2) - value(2,2)*value(3,1);
                    break;
                default:
                    mthrow("Cofact::j_out_of_bound!");
                    break;
            }
            break;
        case 2:
            switch (b) {
                case 1:
                    result = value(1,3)*value(3,2) - value(1,2)*value(3,3);
                    break;
                case 2:
                    result = value(1,1)*value(3,3) - value(1,3)*value(3,1);
                    break;
                case 3:
                    result = value(1,2)*value(3,1) - value(1,1)*value(3,2);
                    break;
                default:
                    mthrow("Cofact::j_out_of_bound!");
                    break;
            }
            break;
        case 3:
            switch (b) {
                case 1:
                    result = value(1,2)*value(2,3) - value(1,3)*value(2,2);
                    break;
                case 2:
                    result = value(1,3)*value(2,1) - value(1,1)*value(2,3);
                    break;
                case 3:
                    result = value(1,1)*value(2,2) - value(1,2)*value(2,1);
                    break;
                default:
                    mthrow("Cofact::j_out_of_bound!");
                    break;
            }
            break;
        default:
            mthrow("Cofact::i_out_of_bound!");
            break;
    }
    return std::move(result);
}

template <class T>
std::ostream& operator<<(std::ostream &os,const Matrix3<T> &m)
{
    os << m[0] << '\n' << m[1] << '\n' << m[2];
    return os;
}

template <typename T>
typename Matrix3<T>::row_type operator*(const Vector3<T>& v, const Matrix3<T>& m)
{
    return typename Matrix3<T>::row_type(m.col(1)*v,m.col(2)*v,m.col(3)*v);
}

template <typename T> void swap(Matrix3<T>& m1,Matrix3<T>& m2){
    using std::swap;
    swap(m1.Mat,m2.Mat);
}
/* ***********************************************************
 End of Matrix3
 *************************************************************/


/* ***********************************************************
 Vector3
 *************************************************************/
/* ***********************************************************
 End of Vector3
 *************************************************************/


#endif
