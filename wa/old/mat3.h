#ifndef MAT3_H
#define MAT3_H

#include "vec3.h"

inline int epsilon(const int i, const int j, const int k);


template <class T>
class Matrix3
{
    //friend ostream& operator<<(ostream &os,const Matrix3<T> &m);
    //friend istream& operator>>(istream &is,Matrix3<T> &m);
    friend class Vector3<T>;
    
    Vector3<Vector3<T> > Mat;

public:

    //initialize
    Matrix3(const T &x11=0,const T &x12=0,const T &x13=0,
            const T &x21=0,const T &x22=0,const T &x23=0,
            const T &x31=0,const T &x32=0,const T &x33=0)
    {   Mat.x.set(x11,x12,x13);Mat.y.set(x21,x22,x23);Mat.z.set(x31,x32,x33);}
    
    Matrix3(const Vector3<T>& l1, const Vector3<T>& l2, const Vector3<T>& l3)
    {   Mat.x=l1; Mat.y=l2; Mat.z=l3;    }

    Matrix3(const Matrix3<T> &obj)
    {   Mat = obj.Mat; }
    
    //Read data
    Vector3<T> row (const int& i, const bool mod= false) const
        //return ith row i(1,2,3)
    {   return Mat[i-1];    }
    
    Vector3<T> col (const int& i) const
        //return ith column i(1,2,3)
    {   return Vector3<T>(Mat[0][i-1],Mat[1][i-1],Mat[2][i-1]); }
    
    T value(const int &i, const int &j,const bool mod= false) const
        //return Mat(i,j) i,j(1,2,3)
    {   return Mat[i-1][j-1]; }
    
    Vector3<T>& operator[] (const int& i)
        //return i+1th row i(0,1,2)
    {   return Mat[i];  }
    
    const Vector3<T>& operator[] (const int& i) const
        //return i+1th row i(0,1,2)
    {   return Mat[i];  }
    
    //arithmetic
    Matrix3<T> operator* (const double& num)
    {   Matrix3<T> tmp=*this; for(int i=0;i<3;++i) tmp[i]*=num; return tmp; }
    
    Matrix3<T> operator* (const int& num)
    {   Matrix3<T> tmp=*this; for(int i=0;i<3;++i) tmp[i]*=num; return tmp; }
    
    Matrix3<T> operator/ (const double& num)
    {   Matrix3<T> tmp=*this; for(int i=0;i<3;++i) tmp[i]/=num; return tmp; }
    
    Matrix3<T> operator/ (const int& num)
    {   Matrix3<T> tmp=*this; for(int i=0;i<3;++i) tmp[i]/=num; return tmp; }
    
    Vector3<T> operator* (const Vector3<T>& v)
    {   return Vector3<T>(Mat[0]*v,Mat[1]*v,Mat[2]*v); }
    
    Matrix3<T> operator* (const Matrix3<T>& m2)
    {   return Matrix3<T>(Vector3<T>(Mat[0]*m2.col(1),Mat[0]*m2.col(2),Mat[0]*m2.col(3)),
                         Vector3<T>(Mat[1]*m2.col(1),Mat[1]*m2.col(2),Mat[1]*m2.col(3)),
                         Vector3<T>(Mat[2]*m2.col(1),Mat[2]*m2.col(2),Mat[2]*m2.col(3))); }
    
    T det() const
    {
        T result=0;
        for(int i=1;i<=3;++i)
            for(int j=1;j<=3;++j)
                for(int k=1;k<=3;k++)
                    result += epsilon(i,j,k)*value(1,i)*value(2,j)*value(3,k);
        return result;
    }
    
    T cofact(const int &i, const int &j) const
        //return Cij i,j(1,2,3)
    {
        vector<T> d4;
        for(int x=1; x<=3; x++)
            for(int y=1; y<=3; y++)
                if( x != i and y != j)
                    d4.push_back(value(x,y));
        assert(d4.size()==4);
        return pow(-1,i+j)*(d4[0]*d4[3]-d4[1]*d4[2]);
    }
    
    Matrix3<T> transposition() const
    {   return Matrix3<T>(col(1),col(2),col(3));    }
    
    Matrix3<T> inverse() const
    {
        assert( det()!=0 );
        Matrix3<T> result(cofact(1,1),cofact(1,2),cofact(1,3),
                          cofact(2,1),cofact(2,2),cofact(2,3),
                          cofact(3,1),cofact(3,2),cofact(3,3));
        return result.transposition()/det();
    }
    
    
    void operator= (const Matrix3<T>& m2)
    {   Mat.x=m2.x; Mat.y=m2.y; Mat.z =m2.z; }
};


template <class T>
istream& operator>>(istream &is,Matrix3<T> &m)
{
    is >> m[0] >> m[1] >> m[2];
    return is;
}

template <class T>
ostream& operator<<(ostream &os,const Matrix3<T> &m)
{
    os << m[0] << '\n' << m[1] << '\n' << m[2];
    return os;
}

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
/*
 T value(const int &i, const int &j,const bool mod= false) const
 //return Mat(i,j) i,j(1,2,3)
 {
 if(mod)
 switch(j%3){
 case 1: return row(i,mod).x;
 case 2: return row(i,mod).y;
 case 0: return row(i,mod).z;}
 else
 switch(j){
 case 1: return row(i,mod).x;
 case 2: return row(i,mod).y;
 case 3: return row(i,mod).z;}
 cerr << "beyond boundary of Matix3!";
 return 0;
 }
 
 Matrix3<T> transposition() const
 {
 Matrix3<T> result=*this;
 for(int x=1; x<=3; x++)
 for(int y=x+1; y<=3; y++)
 {
 result[x-1][y-1] = value(y,x);
 result[y-1][x-1] = value(x,y);
 }
 return result;
 }
 */

#endif


