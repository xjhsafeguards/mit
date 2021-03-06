#ifndef Math_distributionfunction_H
#define Math_distributionfunction_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "Math_const.h"

class Distributionfunction{
    
    //range to calculate distribution function
    double lower,upper;
    int steps;
    double step_length;
    long long data_count = 0;
    int dimension = 1;
    double normalize = 1;
    // add gausian broad parameters
    double _g_sigma;  // as the same unit of step
    
    std::vector<long long> yresult;
    
public:
    Distributionfunction(double in_lower=0, double in_upper=1, int in_steps=100, double in_sig=0): lower(in_lower), upper(in_upper), steps(in_steps),_g_sigma(in_sig){
        set_step_length();
    }
    ~Distributionfunction(){};
    
    //set limit
    inline void set_limit(double in_lower, double in_upper);
    inline void set_upper(double in_upper);
    inline void set_lower(double in_lower);
    //set steps
    inline void set_steps(int in_steps);
    //set dimension of x represent for
    inline void set_dimension(int in_D);
    //set normalized coefficient for y
    inline void set_normalize(double in_norm);
    //set gausian broadening
    inline void set_sigma(double in_sig);
    inline void set_step_sigma(double in_sig);
    
    //read in data to calculate
    void read(const std::vector<double> & data_in);
    void read(typename std::vector<double>::const_iterator begin,typename std::vector<double>::const_iterator end);
    void read(double data_point);
    
    //reset data
    inline void reset();
    inline void reset(double in_lower,double in_upper, int in_steps);
    
    //get result
    std::vector<double> get_x() const;
    std::vector<double> get_y() const;
    const std::vector<long long>& get_ycount() const;
    long long get_valid_count() const;
    
private:
    void set_step_length(){
        assert(steps != 0);
        assert(upper>lower);
        step_length = (upper-lower)/steps;
        reset_yresult();
    }
    void reset_yresult(){
        data_count=0;
        yresult.clear();
        yresult.resize(steps,0);
    }
    void add_value(const double& value){
        if(value >= lower and value < upper)
            ++(yresult[static_cast<int>((value-lower)/step_length)]);
        ++data_count;
    }
    void D1_y(std::vector<double>& result) const;
    void D3_y(std::vector<double>& result) const;
    std::vector<double> g_broaden(const std::vector<double>& in) const;
};

class Distributionfunction2D{
    
    //range to calculate distribution function
    double lower1,upper1,lower2,upper2;
    int steps1,steps2;
    double step_length1,step_length2;
    long long data_count = 0;
    int dimension1 = 1,dimension2 = 1;
    double normalize1 = 1,normalize2 = 1;
    
    std::vector<std::vector<long long>> yresult;
    
public:
    Distributionfunction2D(double in_lower1=0, double in_upper1=1, int in_steps1=100,double in_lower2=0, double in_upper2=1, int in_steps2=100): lower1(in_lower1), upper1(in_upper1), steps1(in_steps1),lower2(in_lower2), upper2(in_upper2), steps2(in_steps2){
        set_step_length();
    }
    
    //set limit
    inline void set_limit1(double in_lower, double in_upper);
    inline void set_upper1(double in_upper);
    inline void set_lower1(double in_lower);
    inline void set_limit2(double in_lower, double in_upper);
    inline void set_upper2(double in_upper);
    inline void set_lower2(double in_lower);
    //set steps
    inline void set_steps1(int in_steps);
    inline void set_steps2(int in_steps);
    //set dimension of x represent for
    inline void set_dimension1(int in_D);
    inline void set_dimension2(int in_D);
    //set normalized coefficient for y
    inline void set_normalize1(double in_norm);
    inline void set_normalize2(double in_norm);
    
    //read in data to calculate
    inline void read(double data1,double data2);
    
    //reset data
    inline void reset();
    inline void reset(double in_lower1,double in_upper1, int in_steps1,double in_lower2,double in_upper2, int in_steps2);
    
    //get result
    std::vector<double> get_x1() const;
    std::vector<double> get_x2() const;
    std::vector<std::vector<double> > get_y() const;
    const std::vector<std::vector<long long> >& get_ycount() const;
    long long get_valid_count() const;
    
private:
    void set_step_length(){
        assert(steps1 != 0);
        assert(steps2 != 0);
        assert(upper1>lower1);
        assert(upper2>lower2);
        step_length1 = (upper1-lower1)/steps1;
        step_length2 = (upper2-lower2)/steps2;
        reset_yresult();
    }
    void reset_yresult(){
        data_count=0;
        yresult.clear();
        yresult.resize(steps1,std::vector<long long>(steps2));
    }
    void add_value(const double& value1, const double& value2){
        if(value1 >= lower1 and value1 < upper1 and value2 >= lower2 and value2 < upper2)
            ++(yresult[static_cast<int>((value1-lower1)/step_length1)][static_cast<int>((value2-lower2)/step_length2)]);
        ++data_count;
    }
    void D11_y(std::vector<std::vector<double> >& result) const;
};

class Distributionfunction3D{
    
    //range to calculate distribution function
    double lower1,upper1,lower2,upper2,lower3,upper3;
    int steps1,steps2,steps3;
    double step_length1,step_length2,step_length3;
    long long data_count = 0;
    //int dimension1 = 1,dimension2 = 1;
    double normalize1 = 1,normalize2 = 1,normalize3 = 1;
    
    std::vector<std::vector<std::vector<long long> > > yresult;
    
public:
    Distributionfunction3D(double in_lower1, double in_upper1, int in_steps1,double in_lower2, double in_upper2, int in_steps2, double in_lower3, double in_upper3, int in_steps3): lower1(in_lower1), upper1(in_upper1), steps1(in_steps1),lower2(in_lower2), upper2(in_upper2), steps2(in_steps2),lower3(in_lower3), upper3(in_upper3), steps3(in_steps3){
        set_step_length();
    }
    Distributionfunction3D(double in_lower=0, double in_upper=1, int in_steps=100): Distributionfunction3D(in_lower,in_upper,in_steps,in_lower,in_upper,in_steps,in_lower,in_upper,in_steps){
        set_step_length();
    }
    
    //set limit
    inline void set_limit1(double in_lower, double in_upper);
    inline void set_upper1(double in_upper);
    inline void set_lower1(double in_lower);
    inline void set_limit2(double in_lower, double in_upper);
    inline void set_upper2(double in_upper);
    inline void set_lower2(double in_lower);
    inline void set_limit3(double in_lower, double in_upper);
    inline void set_upper3(double in_upper);
    inline void set_lower3(double in_lower);
    inline void set_limit(double in_lower, double in_upper);
    inline void set_upper(double in_upper);
    inline void set_lower(double in_lower);
    //set steps
    inline void set_steps1(int in_steps);
    inline void set_steps2(int in_steps);
    inline void set_steps3(int in_steps);
    inline void set_steps(int in_steps);
    //set dimension of x represent for
    //inline void set_dimension1(int in_D);
    //inline void set_dimension2(int in_D);
    //set normalized coefficient for y
    inline void set_normalize1(double in_norm);
    inline void set_normalize2(double in_norm);
    inline void set_normalize3(double in_norm);
    inline void set_normalize(double in_norm);
    //read in data to calculate
    inline void read(double data1,double data2,double data3);
    
    //reset data
    inline void reset();
    inline void reset(double in_lower,double in_upper, int in_steps);
    inline void reset(double in_lower1,double in_upper1, int in_steps1,double in_lower2,double in_upper2, int in_steps2,double in_lower3, double in_upper3, int in_steps3);
    
    //get result
    std::vector<double> get_x1() const;
    std::vector<double> get_x2() const;
    std::vector<double> get_x3() const;
    std::vector<std::vector<std::vector<double> > > get_y() const;
    const std::vector<std::vector<std::vector<long long> > >& get_ycount() const;
    long long get_valid_count() const;
    
private:
    void set_step_length(){
        assert(steps1 != 0);
        assert(steps2 != 0);
        assert(steps3 != 0);
        assert(upper1>lower1);
        assert(upper2>lower2);
        assert(upper3>lower3);
        step_length1 = (upper1-lower1)/steps1;
        step_length2 = (upper2-lower2)/steps2;
        step_length3 = (upper3-lower3)/steps3;
        reset_yresult();
    }
    void reset_yresult(){
        data_count=0;
        yresult.clear();
        yresult.resize(steps1,std::vector<std::vector<long long> >(steps2,std::vector<long long>(steps3)));
    }
    void add_value(const double& value1, const double& value2, const double& value3){
        if(value1 >= lower1 and value1 < upper1 and value2 >= lower2 and value2 < upper2 and value3 >= lower3 and value3 < upper3)
            ++(yresult[static_cast<int>((value1-lower1)/step_length1)][static_cast<int>((value2-lower2)/step_length2)][static_cast<int>((value3-lower3)/step_length3)]);
        ++data_count;
    }
    void D111_y(std::vector<std::vector<std::vector<double> > >& result) const;
};

inline void Distributionfunction::set_limit(double in_lower, double in_upper){
    lower = in_lower;
    upper = in_upper;
    set_step_length();
}
inline void Distributionfunction::set_upper(double in_upper){
    upper = in_upper;
    set_step_length();
}
inline void Distributionfunction::set_lower(double in_lower){
    lower = in_lower;
    set_step_length();
}
inline void Distributionfunction::set_steps(int in_steps){
    steps = in_steps;
    set_step_length();
}
inline void Distributionfunction::set_dimension(int in_D){
    dimension = in_D;
}
inline void Distributionfunction::set_normalize(double in_norm){
    normalize = in_norm;
}
inline void Distributionfunction::set_sigma(double in_sig){
    _g_sigma=in_sig;
}
inline void Distributionfunction::set_step_sigma(double in_sig){
    _g_sigma=in_sig*step_length;
}
inline void Distributionfunction::reset(){
    reset_yresult();
}
inline void Distributionfunction::reset(double in_lower,double in_upper, int in_steps){
    lower = in_lower;
    upper = in_upper;
    steps = in_steps;
    set_step_length();
}
/* Functions for Distrutionfunction2D*/
inline void Distributionfunction2D::set_limit1(double in_lower, double in_upper){
    lower1 = in_lower;
    upper1 = in_upper;
    set_step_length();
}
inline void Distributionfunction2D::set_limit2(double in_lower, double in_upper){
    lower2 = in_lower;
    upper2 = in_upper;
    set_step_length();
}
inline void Distributionfunction2D::set_upper1(double in_upper){
    upper1 = in_upper;
    set_step_length();
}
inline void Distributionfunction2D::set_upper2(double in_upper){
    upper2 = in_upper;
    set_step_length();
}
inline void Distributionfunction2D::set_lower1(double in_lower){
    lower1 = in_lower;
    set_step_length();
}
inline void Distributionfunction2D::set_lower2(double in_lower){
    lower2 = in_lower;
    set_step_length();
}
inline void Distributionfunction2D::set_steps1(int in_steps){
    steps1 = in_steps;
    set_step_length();
}
inline void Distributionfunction2D::set_steps2(int in_steps){
    steps2 = in_steps;
    set_step_length();
}
inline void Distributionfunction2D::set_dimension1(int in_D){
    dimension1 = in_D;
}
inline void Distributionfunction2D::set_dimension2(int in_D){
    dimension2 = in_D;
}
inline void Distributionfunction2D::set_normalize1(double in_norm){
    normalize1 = in_norm;
}
inline void Distributionfunction2D::set_normalize2(double in_norm){
    normalize2 = in_norm;
}
inline void Distributionfunction2D::read(double data1,double data2){
    add_value(data1,data2);
}
inline void Distributionfunction2D::reset(){
    reset_yresult();
}
inline void Distributionfunction2D::reset(double in_lower1,double in_upper1, int in_steps1,double in_lower2,double in_upper2, int in_steps2){
    lower1 = in_lower1;
    upper1 = in_upper1;
    steps1 = in_steps1;
    lower2 = in_lower2;
    upper2 = in_upper2;
    steps2 = in_steps2;
    set_step_length();
}
/* Functions for Distrutionfunction3D*/
inline void Distributionfunction3D::set_limit1(double in_lower, double in_upper){
    lower1 = in_lower;
    upper1 = in_upper;
    set_step_length();
}
inline void Distributionfunction3D::set_limit2(double in_lower, double in_upper){
    lower2 = in_lower;
    upper2 = in_upper;
    set_step_length();
}
inline void Distributionfunction3D::set_limit3(double in_lower, double in_upper){
    lower3 = in_lower;
    upper3 = in_upper;
    set_step_length();
}
inline void Distributionfunction3D::set_limit(double in_lower, double in_upper){
    lower1 = in_lower;
    upper1 = in_upper;
    lower2 = in_lower;
    upper2 = in_upper;
    lower3 = in_lower;
    upper3 = in_upper;
    set_step_length();
}
inline void Distributionfunction3D::set_upper1(double in_upper){
    upper1 = in_upper;
    set_step_length();
}
inline void Distributionfunction3D::set_upper2(double in_upper){
    upper2 = in_upper;
    set_step_length();
}
inline void Distributionfunction3D::set_upper3(double in_upper){
    upper3 = in_upper;
    set_step_length();
}
inline void Distributionfunction3D::set_upper(double in_upper){
    upper1 = in_upper;
    upper2 = in_upper;
    upper3 = in_upper;
    set_step_length();
}
inline void Distributionfunction3D::set_lower1(double in_lower){
    lower1 = in_lower;
    set_step_length();
}
inline void Distributionfunction3D::set_lower2(double in_lower){
    lower2 = in_lower;
    set_step_length();
}
inline void Distributionfunction3D::set_lower3(double in_lower){
    lower3 = in_lower;
    set_step_length();
}
inline void Distributionfunction3D::set_lower(double in_lower){
    lower1 = in_lower;
    lower2 = in_lower;
    lower3 = in_lower;
    set_step_length();
}
inline void Distributionfunction3D::set_steps1(int in_steps){
    steps1 = in_steps;
    set_step_length();
}
inline void Distributionfunction3D::set_steps2(int in_steps){
    steps2 = in_steps;
    set_step_length();
}
inline void Distributionfunction3D::set_steps3(int in_steps){
    steps3 = in_steps;
    set_step_length();
}
inline void Distributionfunction3D::set_steps(int in_steps){
    steps1 = in_steps;
    steps2 = in_steps;
    steps3 = in_steps;
    set_step_length();
}
inline void Distributionfunction3D::set_normalize1(double in_norm){
    normalize1 = in_norm;
}
inline void Distributionfunction3D::set_normalize2(double in_norm){
    normalize2 = in_norm;
}
inline void Distributionfunction3D::set_normalize3(double in_norm){
    normalize3 = in_norm;
}
inline void Distributionfunction3D::set_normalize(double in_norm){
    normalize1 = in_norm;
    normalize2 = in_norm;
    normalize3 = in_norm;
}
inline void Distributionfunction3D::read(double data1,double data2,double data3){
    add_value(data1,data2,data3);
}
inline void Distributionfunction3D::reset(){
    reset_yresult();
}
inline void Distributionfunction3D::reset(double in_lower,double in_upper, int in_steps){
    lower1 = in_lower;
    upper1 = in_upper;
    steps1 = in_steps;
    lower2 = in_lower;
    upper2 = in_upper;
    steps2 = in_steps;
    lower3 = in_lower;
    upper3 = in_upper;
    steps3 = in_steps;
    set_step_length();
}
inline void Distributionfunction3D::reset(double in_lower1,double in_upper1, int in_steps1,double in_lower2,double in_upper2, int in_steps2,double in_lower3, double in_upper3, int in_steps3){
    lower1 = in_lower1;
    upper1 = in_upper1;
    steps1 = in_steps1;
    lower2 = in_lower2;
    upper2 = in_upper2;
    steps2 = in_steps2;
    lower3 = in_lower3;
    upper3 = in_upper3;
    steps3 = in_steps3;
    set_step_length();
}
#endif
