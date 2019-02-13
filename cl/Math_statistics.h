#ifndef Math_statistics_h
#define Math_statistics_h

#include <vector>
#include <set>

class Statistics{
    std::set<double> data;
    
    double data_sum = 0;
    
public:
    Statistics(){}
    ~Statistics(){}
    
    //Input
    inline void read(const std::vector<double> & data_in);
    inline void read(typename std::vector<double>::const_iterator begin,typename std::vector<double>::const_iterator end);
    inline void read(double data_point);
    //actions
    inline double sum() const;
    inline double average() const;
};

inline void Statistics::read(const std::vector<double> & data_in){
    read(data_in.cbegin(),data_in.cend());
}
inline void Statistics::read(typename std::vector<double>::const_iterator begin,typename std::vector<double>::const_iterator end){
    while(begin!=end){
        data.insert(*begin);
        data_sum += *begin;
        ++begin;
    }
}
inline void Statistics::read(double data_point){
    data.insert(data_point);
    data_sum += data_point;
    
}
inline double Statistics::sum() const{
    return data_sum;
}
inline double Statistics::average() const{
    return data_sum/data.size();
}
#endif
