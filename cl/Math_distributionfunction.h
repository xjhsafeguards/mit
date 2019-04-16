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
    
    std::vector<long long> yresult;
    
public:
    Distributionfunction(double in_lower=0, double in_upper=1, int in_steps=100): lower(in_lower), upper(in_upper), steps(in_steps){
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
};
/*
class Distributionfunction2D{
    
    //range to calculate distribution function
    double lower1,upper1,lower2,upper2;
    int steps1,step2;
    double step_length1,step_lenth2;
    long long data_count = 0;
    int dimension1,dimension2 = 1;
    double normalize = 1;
    
    std::vector<std::vector<long long>> yresult;
    
public:
    Distributionfunction(double in_lower1=0, double in_upper1=1, int in_steps1=100,double in_lower2=0, double in_upper2=1, int in_steps2=100): lower1(in_lower1), upper1(in_upper1), steps1(in_steps1),lower2(in_lower2), upper2(in_upper2), steps2(in_steps2){
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
    void read(double data1,double data2);
    
    //reset data
    inline void reset();
    inline void reset(double in_lower1,double in_upper1, int in_steps1,double in_lower2,double in_upper2, int in_steps2);
    
    //get result
    std::vector<double> get_x1() const;
    std::vector<double> get_x2() const;
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
};
   */
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
inline void Distributionfunction::reset(){
    reset_yresult();
}
inline void Distributionfunction::reset(double in_lower,double in_upper, int in_steps){
    lower = in_lower;
    upper = in_upper;
    steps = in_steps;
    set_step_length();
}
#endif
