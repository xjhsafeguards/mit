#ifndef Math_distributionfunction_H
#define Math_distributionfunction_H

#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>

class Distributionfunction{
    
    //range to calculate distribution function
    double lower,upper;
    int steps;
    double step_length;
    int data_count = 0;
    int dimension = 1;
    double normalize = 1;
    
    std::vector<int> yresult;
    
public:
    Distributionfunction(double in_lower=0, double in_upper=1, int in_steps=100): lower(in_lower), upper(in_upper), steps(in_steps){
        set_step_length();
    }
    ~Distributionfunction();
    
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
    void read(typename std::vector<double>::iterator begin,typename std::vector<double>::iterator end);
    
    //get result
    std::vector<double> get_x() const;
    std::vector<double> get_y() const;
    
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
#endif
