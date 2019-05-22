#include "Math_statistics.h"

double Stat_pool::variance() const{
    assert(data.size()>1);
    double tmp_avg = data_sum/data.size();
    double result = 0;
    for( const auto& num : data)
        result += (num-tmp_avg)*(num-tmp_avg);
    return result/(data.size()-1);
}

