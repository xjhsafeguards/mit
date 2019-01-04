#ifndef Math_distributionfunction_H
#define Math_distributionfunction_H

#include <vector>
#include <cassert>

class Distributionfunction{
    
    //range to calculate distribution function
    double lower,upper;
    int steps；
    double step_length;
    std::vector<double> result;
    
public:
    Distributionfunction(double in_lower=0, double in_upper=1, int in_steps=100): lower(in_lower), upper(in_upper), steps(in_steps){
        set_step_length()
    }
    ~Distributionfunction();
    
    //set limit and steps
    void set_upper(double in_upper){
        upper = in_upper;
        set_step_length();
    }
    void set_lower(double in_lower){
        lower = in_lower;
        set_step_length();
    }
    void set_steps(int in_steps){
        steps = in_steps;
        set_step_length();
    }
    
    //read in data to calculate
    void operator()(const std::vector<double> & data_in);
    void operator()(typename std::vector<double>::iterator& begin,typename std::vector<double>::iterator& it end);
    
private:
    void set_step_length(){
        assert(steps != 0);
        step_length = (upper-lower)/steps;
    }
    
}

#endif
