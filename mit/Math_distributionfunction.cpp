#include "Math_distributionfunction.h"

void Distributionfunction::read(const std::vector<double> & data_in){
    read(data_in.cbegin(),data_in.cend());
}
void Distributionfunction::read(typename std::vector<double>::const_iterator begin,typename std::vector<double>::const_iterator end){
    for(;begin!=end;++begin){
        add_value(*begin);
    }
}
std::vector<double> Distributionfunction::get_x() const{
    std::vector<double> result;
    double tmp=lower+step_length/2;
    for(int c=0;c<=steps;++c,tmp+=step_length){
        result.push_back(tmp);
    }
    return std::move(result);
}
std::vector<double> Distributionfunction::get_y() const{
    std::vector<double> result;
    switch(dimension){
        case 1:
            D1_y(result);
            break;
        case 3:
            D3_y(result);
            break;
        default:
            std::cerr << "Distributionfunction_Wrong_Dimension: " << dimension << std::endl;
    }
    return std::move(result);
}
void Distributionfunction::D1_y(std::vector<double>& result) const{
    for(auto it=yresult.cbegin(); it!=yresult.cend(); ++it){
        result.push_back((*it)*normalize/data_count/step_length);
    }
}
void Distributionfunction::D3_y(std::vector<double>& result) const{
    double r=lower;
    for(auto it=yresult.cbegin(); it!=yresult.cend(); ++it,r+=step_length){
        result.push_back((*it)*normalize/data_count/(4*PI*(pow(r+step_length,3)-pow(r,3))/3));
    }
}
