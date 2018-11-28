#ifndef MATH_H
#define MATH_H

#include "gfun.h"

extern const double PI;

//add 2018.4.11 by jianhang
template<typename iterator, typename T>
double fourier_transform(iterator begin, iterator end, T dx, T k, T x0=0)
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
double fourier_transform_function(iterator1 begin, iterator1 end, iterator2 result, double dx, int Ndx, double kmax = -1, double x0=0)
// Ndx hold the number of dx, can be greater than end - begin, when it indicates the term after end is 0
{
    if(end <= begin) throw(std::range_error("FOURIER_TRANSFORM_FUNCTION: begin + 1 >= end"));
    if(Ndx < 0) throw(std::runtime_error("FOURIER_TRANSFORM_FUNCTION: Ndx < 0"));
    double dk = 2*PI/(dx*Ndx);
    int Nkmax = kmax/dk;
    if(kmax < 0) Nkmax = -1; // remove the Nkmax constrain
    for(int i=0; i!= Ndx+1 and i!= Nkmax; i++)
        *(result++) = fourier_transform(begin,end,dx,i*dk,x0);
    return dk;
}

/* derivation */
template<typename iterator1, typename iterator2,typename T>
void derivate(iterator1 begin, iterator1 end, iterator2 result, T dt)
{
    if(end <= begin) throw(std::range_error("VELOCITY: begin >= end"));
    while(begin + 2 < end)
    {
        *(result++) = (*(begin+2) - *begin)/dt/2;
        ++begin;
    }
}


//add 2018.4.5 by jianhang
/* integration */
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

template<typename iterator>
double correlate(iterator begin, iterator end, unsigned step = 0)
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
double correlate_f(iterator begin, iterator end, F f , unsigned step = 0)
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
int correlate_function(iterator1 begin, iterator1 end, iterator2 result, unsigned n=1, unsigned n_jump =1)
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
int correlate_f_function(iterator1 begin, iterator1 end, iterator2 result, F f, unsigned n=1, unsigned n_jump =1)
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

/* Legendre polynomial */
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


//  old ones
inline double Legendre2(double x)
//second order Legendre polynomial
{
    return 0.5*(3*x*x - 1);
}

/* Intergration */

// calculation methods
double Integration(double *y,double *x, int n2, int n1=0);
double Integration(double *y,double dx, int n2, int n1=0);
double Integration(vector<double>& y,double dx, int n2, int n1=0);
double Integration(double *y,double dx, double n2, int n1=0);
//intergration from n1 to n2(index) of list, double means linear interpolation
//n2 > n1

/* Time correlation function */
double Time_correlation(double *y, int n_max, int n_step);
double Time_correlation(const vector<double>&y, int n_max, int n_step);
double Time_correlation(Vector3<double> *y, int n_max, int n_step);
double Time_correlation(const vector<Vector3<double> >& y, int n_max, int n_step);
double Time_correlation(Vector3<double> *y1,Vector3<double> *y2, int n_max, int n_step);
double Time_correlation(Vector3<double> *y1,Vector3<double> *y2,bool *delta, int n_max, int n_step);
//n_max is the last index of list
//n_step is calculate i,i+n_step
//tao = dt*n_step;
double Time_correlation_Legendre2(Vector3<double> *y, int n_max, int n_step);
double Time_correlation_Legendre2(const vector<Vector3<double> >& y, int n_max, int n_step);

double Time_correlation_Legendre(const vector<double>& y, int n_max, int n_step, int n=1);
double Time_correlation_Legendre(const vector<Vector3<double> >& y, int n_max, int n_step, int n=1);
// for n order Legendre

/* Fourier_transform */
double Fourier_transform(double k,double *y,double dx,int n2, int n1=0,double x1=0);
//k contains 2PI

#endif
