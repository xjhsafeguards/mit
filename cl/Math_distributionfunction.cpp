#include "Math_distributionfunction.h"

void Distributionfunction::read(const std::vector<double> & data_in){
    read(data_in.cbegin(),data_in.cend());
}
void Distributionfunction::read(typename std::vector<double>::const_iterator begin,typename std::vector<double>::const_iterator end){
    for(;begin!=end;++begin){
        add_value(*begin);
    }
}
void Distributionfunction::read(double data_point){
    add_value(data_point);
}
std::vector<double> Distributionfunction::get_x() const{
    std::vector<double> result;
    double tmp=lower+step_length/2;
    for(int c=0;c<steps;++c,tmp+=step_length){
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
const std::vector<long long>& Distributionfunction::get_ycount() const{
    return yresult;
}
long long Distributionfunction::get_valid_count() const{
    long long total=0;
    for(const auto& num: yresult)
        total+=num;
    return total;
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

std::vector<double> Distributionfunction2D::get_x1() const{
    std::vector<double> result;
    double tmp=lower1+step_length1/2;
    for(int c=0;c<steps1;++c,tmp+=step_length1){
        result.push_back(tmp);
    }
    return std::move(result);
}
std::vector<double> Distributionfunction2D::get_x2() const{
    std::vector<double> result;
    double tmp=lower2+step_length2/2;
    for(int c=0;c<steps2;++c,tmp+=step_length2){
        result.push_back(tmp);
    }
    return std::move(result);
}
std::vector<std::vector<double> > Distributionfunction2D::get_y() const{
    std::vector<std::vector<double> > result;
    if(dimension1==1 and dimension2==1)
        D11_y(result);
    else
        std::cerr << "Distributionfunction2D_Wrong_Dimension: " << dimension1 << dimension2 << std::endl;
    return std::move(result);
}
const std::vector<std::vector<long long> >& Distributionfunction2D::get_ycount() const{
    return yresult;
}
long long Distributionfunction2D::get_valid_count() const{
    long long total=0;
    for(const auto& list: yresult)
        for(const auto& num: list)
            total+=num;
    return total;
}

void Distributionfunction2D::D11_y(std::vector<std::vector<double> >& result) const{
    for(auto it1=yresult.cbegin(); it1!=yresult.cend(); ++it1){
        std::vector<double> tmp;
        for(auto it2=it1->cbegin(); it2!=it1->cend(); ++it2){
            tmp.push_back((*it2)*normalize1*normalize2/data_count/step_length1/step_length2);
        }
        result.push_back(tmp);
    }
}


std::vector<double> Distributionfunction3D::get_x1() const{
    std::vector<double> result;
    double tmp=lower1+step_length1/2;
    for(int c=0;c<steps1;++c,tmp+=step_length1){
        result.push_back(tmp);
    }
    return result;
}
std::vector<double> Distributionfunction3D::get_x2() const{
    std::vector<double> result;
    double tmp=lower2+step_length2/2;
    for(int c=0;c<steps2;++c,tmp+=step_length2){
        result.push_back(tmp);
    }
    return result;
}
std::vector<double> Distributionfunction3D::get_x3() const{
    std::vector<double> result;
    double tmp=lower3+step_length3/2;
    for(int c=0;c<steps3;++c,tmp+=step_length3){
        result.push_back(tmp);
    }
    return result;
}
std::vector<std::vector<std::vector<double> > > Distributionfunction3D::get_y() const{
    std::vector<std::vector<std::vector<double> > > result;
    D111_y(result);
    return result;
}
const std::vector<std::vector<std::vector<long long> > >& Distributionfunction3D::get_ycount() const{
    return yresult;
}
long long Distributionfunction3D::get_valid_count() const{
    long long total=0;
    for(const auto& list1: yresult)
        for(const auto& list2: list1)
            for(const auto& num: list2)
                total+=num;
    return total;
}

void Distributionfunction3D::D111_y(std::vector<std::vector<std::vector<double> > >& result) const{
    for(auto it1=yresult.cbegin(); it1!=yresult.cend(); ++it1){
        std::vector<std::vector<double> > tmp1;
        for(auto it2=it1->cbegin(); it2!=it1->cend(); ++it2){
            std::vector<double> tmp2;
            for(auto it3=it2->cbegin(); it3!=it2->cend(); ++it2)
                tmp2.push_back((*it3)*normalize1*normalize2*normalize3/data_count/step_length1/step_length2/step_length3);
            tmp1.push_back(tmp2);
        }
        result.push_back(tmp1);
    }
}
