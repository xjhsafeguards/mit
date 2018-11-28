#include "math.h"

const double PI = 3.141592653589793238460;

double Integration(double *y,double *x, int n2, int n1)
{
    assert(n1 >= 0);
    assert(n2 > n1);
    
    double result = 0;
    result += 0.5*y[n1]*(x[n1+1]-x[n1]);
    cout << result << endl;
    for(int i=n1+1; i<n2; i++)
    {
        result += 0.5*y[i]*(x[i+1]-x[i-1]);
        cout << result << endl;
    }
    result += 0.5*y[n2]*(x[n2]-x[n2-1]);
    cout << result << endl;
    
    return result;
}

double Integration(double *y,double dx, int n2, int n1)
{
    assert(n2 > n1);
    
    double result = 0;
    result += 0.5*y[n1];
    for(int i=n1+1; i<n2; i++)
        result += y[i];
    result += 0.5*y[n2];
    result *= dx;
    
    return result;
}

double Integration(vector<double>& y,double dx, int n2, int n1)
{
    assert(n2 > n1);
    
    double result = 0;
    result += 0.5*y[n1];
    for(int i=n1+1; i<n2; i++)
        result += y[i];
    result += 0.5*y[n2];
    result *= dx;
    
    return result;
}

double Integration(double *y,double dx, double n2, int n1)
{
    int n = (int)n2;
    double d = n2 - n;
    return Integration(y,dx,n,n1)+d*dx*(y[n]+(y[n+1]-y[n])*d/2);
}

double Time_correlation(double *y, int n_max, int n_step)
{
    assert(n_max>n_step);
    assert(n_step>=0);
    
    double *y1;
    y1 = new double [n_max-n_step+1];
    for(int i=0; i<(n_max-n_step+1); i++)
        y1[i]=y[i]*y[i+n_step];
    
    double result = Integration(y1,1.0/(n_max-n_step),n_max-n_step);
    delete [] y1;
    return result;
    
}

double Time_correlation(const vector<double>& y, int n_max, int n_step)
{
    assert(n_max>n_step);
    assert(n_step>=0);
    
    double *y1;
    y1 = new double [n_max-n_step+1];
    for(int i=0; i<(n_max-n_step+1); i++)
        y1[i]=y[i]*y[i+n_step];
    
    double result = Integration(y1,1.0/(n_max-n_step),n_max-n_step);
    delete [] y1;
    return result;
    
}

double Time_correlation(Vector3<double> *y, int n_max, int n_step)
{
    assert(n_max>n_step);
    assert(n_step>=0);
    
    double *y1;
    y1 = new double [n_max-n_step+1];
    for(int i=0; i<(n_max-n_step+1); i++)
        y1[i]=y[i]*y[i+n_step];
    double result = Integration(y1,1.0/(n_max-n_step),n_max-n_step);
    delete [] y1;
    return result;
}

double Time_correlation(const vector<Vector3<double> >& y, int n_max, int n_step)
{
    assert(n_max>n_step);
    assert(n_step>=0);
    
    double *y1;
    y1 = new double [n_max-n_step+1];
    for(int i=0; i<(n_max-n_step+1); i++)
        y1[i]=y[i]*y[i+n_step];
    double result = Integration(y1,1.0/(n_max-n_step),n_max-n_step);
    delete [] y1;
    return result;
}

double Time_correlation(Vector3<double> *y1,Vector3<double> *y2, int n_max, int n_step)
{
    assert(n_max>n_step);
    assert(n_step>=0);
    
    double *y0;
    y0 = new double [n_max-n_step+1];
    for(int i=0; i<(n_max-n_step+1); i++)
        y0[i]=y1[i]*y2[i+n_step];
    double result = Integration(y0,1.0/(n_max-n_step),n_max-n_step);
    delete [] y0;
    return result;
}

double Time_correlation(Vector3<double> *y1,Vector3<double> *y2,bool *delta, int n_max, int n_step)
{
    assert(n_max>n_step);
    assert(n_step>=0);
    
    double *y0;
    y0 = new double [n_max-n_step+1];
    for(int i=0; i<(n_max-n_step+1); i++)
    {
        if(delta[i])
            y0[i]=y1[i]*y2[i+n_step];
        else
            y0[i]=0;
    }
    double result = Integration(y0,1.0/(n_max-n_step),n_max-n_step);
    delete [] y0;
    return result;
}

double Time_correlation_Legendre2(Vector3<double> *y, int n_max, int n_step)
{
    assert(n_max>n_step);
    assert(n_step>=0);
    
    double *y1;
    y1 = new double [n_max-n_step+1];
    for(int i=0; i<(n_max-n_step+1); i++)
        y1[i]= Legendre2(y[i]*y[i+n_step]);
    double result=Integration(y1,1.0/(n_max-n_step),n_max-n_step);
    delete [] y1;
    return result;
}

double Time_correlation_Legendre2(const vector<Vector3<double> >& y, int n_max, int n_step)
{
    assert(n_max>n_step);
    assert(n_step>=0);
    
    double *y1;
    y1 = new double [n_max-n_step+1];
    for(int i=0; i<(n_max-n_step+1); i++)
        y1[i]= Legendre2(y[i]*y[i+n_step]);
    double result=Integration(y1,1.0/(n_max-n_step),n_max-n_step);
    delete [] y1;
    return result;
}

double Time_correlation_Legendre(const vector<double>& y, int n_max, int n_step, int n)
{
    assert(n_max>n_step);
    assert(n_step>=0);
    
    double *y1;
    y1 = new double [n_max-n_step+1];
    for(int i=0; i<(n_max-n_step+1); i++)
        y1[i]= Legendre(n,y[i]*y[i+n_step]);
    double result=Integration(y1,1.0/(n_max-n_step),n_max-n_step);
    delete [] y1;
    return result;
}

double Time_correlation_Legendre(const vector<Vector3<double> >& y, int n_max, int n_step, int n)
{
    assert(n_max>n_step);
    assert(n_step>=0);
    
    double *y1;
    y1 = new double [n_max-n_step+1];
    for(int i=0; i<(n_max-n_step+1); i++)
        y1[i]= Legendre(n,y[i]*y[i+n_step]);
    double result=Integration(y1,1.0/(n_max-n_step),n_max-n_step);
    delete [] y1;
    return result;
}

double Fourier_transform(double k,double *y,double dx,int n2, int n1,double x1)
{
    assert(n2 > n1);
    
    double result = 0,tmp;
    
    double *y1;
    y1 = new double [n2+1];
    for(int i=0; i<(n2+1); i++)
        y1[i]=y[i]*cos(k*(dx*i+x1));
    
    result = Integration(y1,dx,n2,n1);
    delete [] y1;
    return result;
    
    /*
    tmp = Integration(y1,dx,n2,n1);
    result += tmp*tmp;
    
    for(int i=0; i<(n2+1); i++)
        y1[i]=y[i]*sin(k*(dx*i+x1));
    tmp = Integration(y1,dx,n2,n1);
    result += tmp*tmp;
    
    return sqrt(result);
     */
}

