#ifndef MFUNC_H
#define MFUNC_H


namespace math_fucntion {
    
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

}

using math_fucntion::Legendre;
using math_fucntion::Legendre2;
using math_fucntion::epsilon;

#endif
