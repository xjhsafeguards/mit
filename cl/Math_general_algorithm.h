#ifndef MATH_GENERAL_ALGORITHM_H
#define MATH_GENERAL_ALGORITHM_H

//#include <cassert>
#include <stdexcept>

namespace MGA{
    
    // sort
    
    // Fx stores in seq container begin to end
    template<typename iterator, typename T> double fourier_transform(iterator begin, iterator end, T dx, T k, T x0=0);
    template<typename iterator1, typename iterator2> double fourier_transform_function(iterator1 begin, iterator1 end, iterator2 result, double dx, int Ndx, double kmax = -1, double x0=0);
    template<typename iterator1, typename iterator2,typename T> void derivate(iterator1 begin, iterator1 end, iterator2 result, T dt);
    template<typename iterator, typename T> double integrate(iterator begin, iterator end, T dx);
    template<typename iterator, typename T> double integrate_simpson(iterator begin, iterator end, T dx);
    
    // Fx stores in F(x)
    
    
    
    /************************************************************
     fourier_transform
     *************************************************************/
    // f(x0+dx) stores in {begin to end}, calculate Fourier of f at point k
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
    // f(x0+dx) stores in {begin to end}, result could be a back inserter, calculate Fourier of f at point 0+n*dk
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
    
    // f(x0+dx) stores in {begin to end}
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
    
    
    /******************************************************/
    
    //integrate F(x) from x0 to x0+N*dx with step dx
    template<typename F>
    double function_integrate(double x0, double dx, int N, F in_f)
    {
        if(N<2) throw(std::range_error("INTEGRATE: lest than two point"));
        double result= -0.5*F(x0);
        for(int i=0;i<N;++i){
            result += in_f(x0);
            x0 += dx;
        }
        result -=  -0.5*in_f(x0);
        return result*dx;
    }
    
    //integrate F(x,y) from x0 to x0+N*dx with step dx and return the result of G(y) in a seq container
    template<typename F,typename iterator>
    void partial_integrate(double x0, double dx, int N, double y0, double dy, int Ny, iterator result, F in_f)
    {
        for(int i=0; i<N; ++i){
            *(result++) = function_integrate(x0,dx,N,[&y0](double x){return F(x,y0);});
            y0 += dy;
        }
    }
}

using namespace MGA;

#endif

