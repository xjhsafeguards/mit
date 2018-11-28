#ifndef GMATH_H
#define GMATH_H

template<typename> class Vector3;
template <typename> class Matrix3;

#include "gfun.h"

extern const double PI;
double Legendre(int, double);
double Legendre2(double);
int epsilon(const int, const int, const int);
template<typename> class Vector3;
template<typename T> istream& operator>>(istream&, Vector3<T>&);
template<typename T> ostream& operator<<(ostream&, const Vector3<T>&);
template <typename> class Matrix3;
template <typename T> istream& operator>>(istream&, Matrix3<T>&);
template <typename T> ostream& operator<<(ostream&, const Matrix3<T>&);
template <typename T> typename Matrix3<T>::row_type operator*(const Vector3<T>&, const Matrix3<T>&);

template<typename iterator, typename T> double fourier_transform(iterator begin, iterator end, T dx, T k, T x0=0);
template<typename iterator1, typename iterator2> double fourier_transform_function(iterator1 begin, iterator1 end, iterator2 result, double dx, int Ndx, double kmax = -1, double x0=0);
template<typename iterator1, typename iterator2,typename T> void derivate(iterator1 begin, iterator1 end, iterator2 result, T dt);
template<typename iterator, typename T> double integrate(iterator begin, iterator end, T dx);
template<typename iterator, typename T> double integrate_simpson(iterator begin, iterator end, T dx);
template<typename iterator> double correlate(iterator begin, iterator end, unsigned step = 0);
template<typename iterator, typename F> double correlate_f(iterator begin, iterator end, F f , unsigned step = 0);
template<typename iterator1, typename iterator2> int correlate_function(iterator1 begin, iterator1 end, iterator2 result, unsigned n=1, unsigned n_jump =1);
template<typename iterator1, typename iterator2, typename F> int correlate_f_function(iterator1 begin, iterator1 end, iterator2 result, F f, unsigned n=1, unsigned n_jump =1);
template<typename iterator> double correlate2(iterator begin1, iterator end1, iterator begin2, iterator end2, unsigned step = 0);
template<typename iterator, typename F> double correlate2_f(iterator begin1, iterator end1, iterator begin2, iterator end2, F f , unsigned step = 0);
template<typename iterator1, typename iterator2> int correlate2_function(iterator1 begin1, iterator1 end1, iterator1 begin2, iterator1 end2, iterator2 result, unsigned n=1, unsigned n_jump =1);
template<typename iterator1, typename iterator2, typename F> int correlate2_f_function(iterator1 begin1, iterator1 end1, iterator1 begin2, iterator1 end2, iterator2 result, F f, unsigned n=1, unsigned n_jump =1);

/* ***********************************************************
 Legendre polynomial
 *************************************************************/

//nth order Legendre polynomial
inline double Legendre(int n, double x)
{
    switch(n)
    {
        case 0:
            return 1.0;
        case 1:
            return x;
        default:
            return ((2*n-1)*x*Legendre(n-1,x) - (n-1)*Legendre(n-2,x))/n;
    }
}

//second order Legendre polynomial
inline double Legendre2(double x)
{
    return 0.5*(3*x*x - 1);
}

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
    friend istream& operator>> <T>(istream&, Vector3<T>&);
    friend ostream& operator<< <T>(ostream&, const Vector3<T>&);
    template<typename T2> friend class Matrix3;
    
    T x;
    T y;
    T z;
    
public:
    //type
    typedef T value_type;
    //constructor
    Vector3(const T& x1=0, const T& x2=0, const T& x3=0): x(x1),y(x2),z(x3){}
    Vector3(const Vector3<T> &obj): x(obj.x),y(obj.y),z(obj.z){}
    //assignment
    void set(const T &x1,const T &x2, const T &x3)
    {
        x = x1;
        y = x2;
        z = x3;
    }
    Vector3& operator= (const Vector3 &v2)
    {
        this->set(v2.x,v2.y,v2.z);
        return *this;
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
    T distance(const Vector3 &v2) const
    {
        return (*this-v2).norm();
    }
    // calculate the angle between this-v2 , return the anlge in degree
    double vangle(const Vector3 &v2) const
    {
        return acos(this->dot(v2)/this->norm()/v2.norm())/PI*180;
    }
    
    //old functions
    
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
        return acos((d12*d12 + d13*d13 - d23*d23)/(2*d12*d13))/PI*180;
    }
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
    
private:
    void vthrow(string information) const
    {
        throw(runtime_error("VECTOR3::" + information));
    }
    
};

template <typename T>
istream& operator>>(istream& is, Vector3<T> &vec)
{
    is >> vec.x >> vec.y >> vec.z;
    return is;
}

template <typename T>
ostream& operator<<(ostream& os, const Vector3<T> &vec)
{
    os << setw(os.precision()+4) << vec.x << " " << setw(os.precision()+4) << vec.y << " " << setw(os.precision()+4) << vec.z;
    return os;
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
    
    friend ostream& operator<< <T>(ostream&, const Matrix3<T>&);
    friend istream& operator>> <T>(istream&, Matrix3<T>&);
    friend Matrix3<T>::row_type operator*<T>(const Vector3<T>&, const Matrix3<T>&);
    
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
    //assignment
    Matrix3& set(const T &x11=0,const T &x12=0,const T &x13=0,
                const T &x21=0,const T &x22=0,const T &x23=0,
                const T &x31=0,const T &x32=0,const T &x33=0)
    {
        Mat = Vector3<Vector3<T> >(Vector3<T>(x11,x12,x13),Vector3<T>(x21,x22,x23),Vector3<T>(x31,x32,x33));
        return *this;
    }
    Matrix3& operator= (const Matrix3<T>& m2)
    {
        Mat = m2.Mat;
        return *this;
    }
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
        return tmp;
    }
    Matrix3 operator/ (const double& num) const
    {
        Matrix3<T> tmp=*this;
        for(int i=0;i<3;++i)
            tmp[i]/=num;
        return tmp;
    }
    col_type operator* (const Vector3<T>& v) const
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
        return result;
    }
    T cofact(const int &i, const int &j) const     //return Cij i,j(1,2,3)
    {
        vector<T> d4;
        for(int x=1; x<=3; x++)
            for(int y=1; y<=3; y++)
                if( x != i and y != j)
                    d4.push_back(value(x,y));
        assert(d4.size()==4);
        return pow(-1,i+j)*(d4[0]*d4[3]-d4[1]*d4[2]);
    }
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

template <typename T>
typename Matrix3<T>::row_type operator*(const Vector3<T>& v, const Matrix3<T>& m)
{
    return typename Matrix3<T>::row_type(m.col(1)*v,m.col(2)*v,m.col(3)*v);
}
/* ***********************************************************
 End of Matrix3
 *************************************************************/

/* ***********************************************************
 fourier_transform
 *************************************************************/
template<typename iterator, typename T>
double fourier_transform(iterator begin, iterator end, T dx, T k, T x0)
{
    if(end <= begin) throw(std::range_error("FOURIER_TRANSFORM: begin >= end"));
    double result = -0.5*(*begin)*cos(k*x0);
    while(begin != end)
    {
        result += *(begin++) * cos(k*x0);
        x0 += dx;
    }
    result -= 0.5* (*(end-1)) * cos(k*(x0-dx));
    return result*dx;
}

template<typename iterator1, typename iterator2>
double fourier_transform_function(iterator1 begin, iterator1 end, iterator2 result, double dx, int Ndx, double kmax, double x0)
// Ndx hold the number of dx, can be greater than end - begin, when it indicates the term after end is 0
{
    if(end <= begin) throw(std::range_error("FOURIER_TRANSFORM_FUNCTION: begin + 1 >= end"));
    if(Ndx < 0) throw(std::runtime_error("FOURIER_TRANSFORM_FUNCTION: Ndx < 0"));
    double dk = 2*PI/(dx*Ndx);
    int Nkmax = kmax/dk;
    if(kmax < 0) Nkmax = -2; // remove the Nkmax constrain
    for(int i=0; i!= Ndx+1 and i!= Nkmax+1; i++)
        *(result++) = fourier_transform(begin,end,dx,i*dk,x0);
    return dk;
}

/* ***********************************************************
 derivation
 *************************************************************/

template<typename iterator1, typename iterator2,typename T>
void derivate(iterator1 begin, iterator1 end, iterator2 result, T dt)
{
    if(end <= begin) throw(std::range_error("DERIVATION: begin >= end"));
    while(begin + 2 < end)
    {
        *(result++) = (*(begin+2) - *begin)/dt/2;
        ++begin;
    }
}


/* ***********************************************************
 integration
 *************************************************************/

template<typename iterator, typename T>
double integrate(iterator begin, iterator end, T dx)
{
    //typedef typename iterator_traits<iterator>::value_type valuetype;
    if(end <= begin) throw(std::range_error("INTEGRATE: begin >= end"));
    double result = -0.5*(*begin+*(end-1));
    while(begin != end)
        result += *(begin++);
    return result*dx;
}

template<typename iterator, typename T>
double integrate_simpson(iterator begin, iterator end, T dx)
{
    //typedef typename iterator_traits<iterator>::value_type valuetype;
    if(end <= begin) throw(std::range_error("INTEGRATE: begin >= end"));
    double result = 0;
    while(begin != end)
    {
        if(begin + 2 < end)
        {
            result += *begin;
            result += *(++begin)*4;
            result += *(++begin);
        }
        else if(begin + 2 == end)
        {
            result += *begin * 3/2;
            result += *(++begin) * 3/2;
        }
        else
            ++begin;
    }
    return result*dx/3;
}

/* ***********************************************************
 correlate
 *************************************************************/
template<typename iterator>
double correlate(iterator begin, iterator end, unsigned step)
{
    vector<double> temp;
    iterator begin2 = begin + step;
    if(end <= begin2) throw(std::range_error("CORRELATE: begin + step >= end"));
    temp.reserve(end-begin2);
    while(begin2 != end)
        temp.push_back(*(begin++)* *(begin2++));
    return (temp.size()==1)? temp[0] : integrate(temp.cbegin(),temp.cend(),1.0/(temp.size()-1));
}

template<typename iterator, typename F>
double correlate_f(iterator begin, iterator end, F f , unsigned step)
{
    vector<double> temp;
    iterator begin2 = begin + step;
    if(end <= begin2) throw(std::range_error("CORRELATE: begin + step >= end"));
    temp.reserve(end-begin2);
    while(begin2 != end)
        temp.push_back(f(*(begin++), *(begin2++)));
    return (temp.size()==1)? temp[0] : integrate(temp.cbegin(),temp.cend(),1.0/(temp.size()-1));
}

template<typename iterator1, typename iterator2>
int correlate_function(iterator1 begin, iterator1 end, iterator2 result, unsigned n, unsigned n_jump)
// total size of output n, jump between output n_jump
{
    if( n_jump == 0) throw(std::runtime_error("CORRELATE_FUNCTION: n_jump = 0"));
    if( end - begin <= 1) throw(std::runtime_error("CORRELATE_FUNCTION: end - begin <= 1"));
    int i = 0;
    for(; i<n and begin + i*n_jump < end; i++)
    {
        *(result++) = correlate(begin,end,i*n_jump);
    }
    return i;
}

template<typename iterator1, typename iterator2, typename F>
int correlate_f_function(iterator1 begin, iterator1 end, iterator2 result, F f, unsigned n, unsigned n_jump)
// total size of output n, jump between output n_jump
{
    if( n_jump == 0) throw(std::runtime_error("CORRELATE_FUNCTION: n_jump = 0"));
    if( end - begin <= 1) throw(std::runtime_error("CORRELATE_FUNCTION: end - begin <= 1"));
    int i = 0;
    for(; i<n and begin + i*n_jump < end; i++)
    {
        *(result++) = correlate_f(begin,end,f,i*n_jump);
    }
    return i;
}

template<typename iterator>
double correlate2(iterator begin1, iterator end1, iterator begin2, iterator end2, unsigned step)
{
    vector<double> temp;
    begin2 = begin2 + step;
    if(end2 <= begin2) throw(std::range_error("CORRELATE2: begin2 + step >= end"));
    if(end1 <= begin1) throw(std::range_error("CORRELATE2: begin1 >= end1"));
    temp.reserve(end2-begin2);
    while(begin2 != end2 and begin1!=end1)
        temp.push_back(*(begin1++)* *(begin2++));
    return (temp.size()==1)? temp[0] : integrate(temp.cbegin(),temp.cend(),1.0/(temp.size()-1));
}

template<typename iterator, typename F>
double correlate2_f(iterator begin1, iterator end1, iterator begin2, iterator end2, F f , unsigned step)
{
    vector<double> temp;
    begin2 = begin2 + step;
    if(end2 <= begin2) throw(std::range_error("CORRELATE2: begin2 + step >= end"));
    if(end1 <= begin1) throw(std::range_error("CORRELATE2: begin1 >= end1"));
    temp.reserve(end2-begin2);
    while(begin2 != end2 and begin1!=end1)
        temp.push_back(f(*(begin1++), *(begin2++)));
    return (temp.size()==1)? temp[0] : integrate(temp.cbegin(),temp.cend(),1.0/(temp.size()-1));
}

template<typename iterator1, typename iterator2>
int correlate2_function(iterator1 begin1, iterator1 end1, iterator1 begin2, iterator1 end2, iterator2 result, unsigned n, unsigned n_jump)
// total size of output n, jump between output n_jump
{
    if( n_jump == 0) throw(std::runtime_error("CORRELATE_FUNCTION2: n_jump = 0"));
    if( end1 - begin1 <= 1) throw(std::runtime_error("CORRELATE_FUNCTION2: end1 - begin1 <= 1"));
    if( end2 - begin2 <= 1) throw(std::runtime_error("CORRELATE_FUNCTION2: end2 - begin2 <= 1"));
    int i = 0;
    for(; i<n and begin1 + i*n_jump < end1 and begin2 + i*n_jump < end2 ; i++)
    {
        *(result++) = correlate2(begin1,end1,begin2,end2,i*n_jump);
    }
    return i;
}

template<typename iterator1, typename iterator2, typename F>
int correlate2_f_function(iterator1 begin1, iterator1 end1, iterator1 begin2, iterator1 end2, iterator2 result, F f, unsigned n, unsigned n_jump)
// total size of output n, jump between output n_jump
{
    if( n_jump == 0) throw(std::runtime_error("CORRELATE_FUNCTION2: n_jump = 0"));
    if( end1 - begin1 <= 1) throw(std::runtime_error("CORRELATE_FUNCTION2: end1 - begin1 <= 1"));
    if( end2 - begin2 <= 1) throw(std::runtime_error("CORRELATE_FUNCTION2: end2 - begin2 <= 1"));
    int i = 0;
    for(; i<n and begin1 + i*n_jump < end1 and begin2 + i*n_jump < end2 ; i++)
    {
        *(result++) = correlate2_f(begin1,end1,begin2,end2,f,i*n_jump);
    }
    return i;
}

// add 2018.5.21 ny jianhang
template<typename iterator,  typename iteratorc, typename Fc>
double correlate3(iterator begin1, iterator end1, iterator begin2, iterator end2, iteratorc cbegin1, iteratorc cbegin2, Fc fc,  unsigned step)
{
    vector<double> temp;
    begin2 = begin2 + step;
    cbegin2 = cbegin2 + step;
    assert(end2 > begin2);
    assert(end1 > begin1);
    temp.reserve(end2-begin2);
    while(begin2 != end2 and begin1!=end1){
        auto fcv = fc( *(cbegin1++),*(cbegin2++));
        if(fcv)
            temp.push_back(*(begin1)* *(begin2) * fcv) ;
        else
            temp.push_back(0);
        begin1++;
        begin2++;
    }
    return (temp.size()==1)? temp[0] : integrate(temp.cbegin(),temp.cend(),1.0/(temp.size()-1));
}

template<typename iterator1, typename iterator2, typename iteratorc, typename Fc>
int correlate3_function(iterator1 begin1, iterator1 end1, iterator1 begin2, iterator1 end2, iterator2 result,  iteratorc cbegin1, iteratorc cbegin2, Fc fc, unsigned n, unsigned n_jump)
// total size of output n, jump between output n_jump
{
    if( n_jump == 0) throw(std::runtime_error("CORRELATE_FUNCTION2: n_jump = 0"));
    if( end1 - begin1 <= 1) throw(std::runtime_error("CORRELATE_FUNCTION2: end1 - begin1 <= 1"));
    if( end2 - begin2 <= 1) throw(std::runtime_error("CORRELATE_FUNCTION2: end2 - begin2 <= 1"));
    int i = 0;
    for(; i<n and begin1 + i*n_jump < end1 and begin2 + i*n_jump < end2 ; i++)
    {
        *(result++) = correlate3(begin1,end1,begin2,end2,cbegin1,cbegin2,fc,i*n_jump);
    }
    return i;
}
/* ***********************************************************
 Vector3
 *************************************************************/
/* ***********************************************************
 End of Vector3
 *************************************************************/

#endif
