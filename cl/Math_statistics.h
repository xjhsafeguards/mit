#ifndef Math_statistics_h
#define Math_statistics_h

#include <cassert>
#include <algorithm>
#include <vector>
#include <set>
#include <math.h>

class Stat_pool{
    std::multiset<double> data;
    
    double data_sum = 0;
    
public:
    Stat_pool(){}
    ~Stat_pool(){}
    
    //Input
    inline void read(double data_point);
    template <typename InputIterator>
    void read(InputIterator begin,InputIterator end);
    inline void read(const std::vector<double> & data_in);
    inline void read(const Stat_pool& data_in);
    //actions
    inline double count() const;
    inline double sum() const;
    inline double average() const;
    inline double ave() const;
    
    inline double max() const;
    inline double min() const;
    
    double variance() const;
    double standard_deviation() const;
    inline double sd() const;
    
};

template <typename InputIterator>
void Stat_pool::read(InputIterator begin,InputIterator end){
    while(begin!=end){
        data.insert(*begin);
        data_sum += *begin;
        ++begin;
    }
}
inline void Stat_pool::read(double data_point){
    data.insert(data_point);
    data_sum += data_point;
    
}
inline void Stat_pool::read(const std::vector<double> & data_in){
    read(data_in.cbegin(),data_in.cend());
}
inline void Stat_pool::read(const Stat_pool& data_in){
    read(data_in.data.cbegin(),data_in.data.cbegin());
}
inline double Stat_pool::count() const{
    return data.size();
}
inline double Stat_pool::sum() const{
    return data_sum;
}
inline double Stat_pool::average() const{
    assert(data.size()!=0);
    return data_sum/data.size();
}
inline double Stat_pool::ave() const{
    assert(data.size()!=0);
    return data_sum/data.size();
}
inline double Stat_pool::max() const{
    return *std::max_element(data.cbegin(),data.cend());
}
inline double Stat_pool::min() const{
    return *std::min_element(data.cbegin(),data.cend());
}
inline double Stat_pool::standard_deviation() const{
    return sqrt(variance());
}
inline double Stat_pool::sd() const{
    return standard_deviation();
}
#endif
