#ifndef MATH_H
#define MATH_H

#include "Math_const.h"
#include "Math_linearalgebra.h"
#include "Math_distributionfunction.h"
#include "Math_statistics.h"
#include "Math_hungarian.h"
#include "Math_optimal_rotation.h"
#include "Math_general_algorithm.h"

#endif

#ifndef MATH_H
extern const double PI;
double Legendre(int, double);
double Legendre2(double);
int epsilon(const int, const int, const int);
template<typename> class Vector3;
template<typename T> std::istream& operator>>(std::istream&, Vector3<T>&);
template<typename T> std::ostream& operator<<(std::ostream&, const Vector3<T>&);
template <typename> class Matrix3;
template <typename T> std::istream& operator>>(std::istream&, Matrix3<T>&);
template <typename T> std::ostream& operator<<(std::ostream&, const Matrix3<T>&);
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


/************************************************************
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
    std::vector<double> temp;
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
    std::vector<double> temp;
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
    std::vector<double> temp;
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
    std::vector<double> temp;
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
    std::vector<double> temp;
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

#endif
