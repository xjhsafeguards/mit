#include <iomanip>

#include "Math_distributionfunction.h"
using namespace std;

template<typename T>
void Print(const vector<T>& inv,ostream& os=cout){
    for( const auto& d: inv)
        os << setw(10) << d;
    os << endl;
}


int main(){

}

/* Distributionfunction test;
 Distributionfunction DF(0,10,10);
 vector<double> v,v2,v3;
 for(int i=0;i!=10000000;++i)
 v.push_back(10.0*rand()/double(RAND_MAX));
 DF.read(v);
 Print(DF.get_x());
 Print(DF.get_y());
 
 DF.set_upper(100);
 v2 = vector<double>({1,2,3,4,5,6,7});
 DF.read(v2);
 DF.set_normalize(10);
 Print(DF.get_x());
 Print(DF.get_y());
 
 DF.reset(0,1,10);
 Print(DF.get_x());
 Print(DF.get_ycount());
 Print(DF.get_y());
 DF.set_dimension(3);
 for(int i=0;i!=10000000;++i)
 v3.push_back(sqrt(pow(rand()/double(RAND_MAX),2)+pow(rand()/double(RAND_MAX),2)+pow(rand()/double(RAND_MAX),2)));
 DF.read(v3);
 Print(DF.get_x());
 Print(DF.get_ycount());
 Print(DF.get_y());
 */
