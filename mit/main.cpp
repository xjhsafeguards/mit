#include <iomanip>
#include <chrono>

#include "Math.h"
#include "Cell_position.h"

using namespace std;
using namespace chrono;


template<typename T>
void Print(const vector<T>& inv,ostream& os=cout){
    for( const auto& d: inv)
        os << setw(10) << d;
    os << endl;
}


int main(){
    
    position* p1,*p2,*p3;
    box bo(Matrix3<double>(2,0,0,0,3,0,0,0,4));
    auto b=make_shared<box>(bo);
    p1 = new frac_position(Vector3<double>(0.1,0.5,0.2),b);
    p2 = new frac_position(Vector3<double>(0.2,0.5,0.2),b);
    p3 = new frac_position(Vector3<double>(0.1,0.4,0.2),b);
    cout << p1->distance(*p2) << endl;
    auto start = system_clock::now();
    cout << p1->angle(*p2,*p3) << endl;
    auto end   = system_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout <<  "used "
    << double(duration.count()) * microseconds::period::num / microseconds::period::den
    << " second " << endl;
    
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

/* Math_linearalgebra
Vector3<double> v1(1.5,3,8),v2(2,4,7);
Matrix3<double> m1(2,0,0,0,3,0,0,0,4),m2(v1,v2,Vector3<double>(2,5,2));
cout << v1 << endl;
cout << "norm: " << v1.norm() << endl;
cout << "uniform: " << v1.uniform() << endl;
swap(v1,v2);
cout << v1 << endl;
cout << m1 << endl;
cout << m2 << endl;
auto start = system_clock::now();
cout << m2.inverse() << "  ";
auto end   = system_clock::now();
auto duration = duration_cast<microseconds>(end - start);
cout <<  "used "
<< double(duration.count()) * microseconds::period::num / microseconds::period::den
<< " second " << endl;

start = system_clock::now();
cout << m2.inverse2() << "  ";
end   = system_clock::now();
duration = duration_cast<microseconds>(end - start);
cout <<  "used "
<< double(duration.count()) * microseconds::period::num / microseconds::period::den
<< " second " << endl;
*/

/*
 position* p1,*p2,*p3;
 box bo(Matrix3<double>(2,0,0,0,3,0,0,0,4));
 auto b=make_shared<box>(bo);
 p1 = new frac_position(Vector3<double>(0.1,0.5,0.2),b);
 p2 = new frac_position(Vector3<double>(0.2,0.5,0.2),b);
 p3 = new frac_position(Vector3<double>(0.1,0.4,0.2),b);
 cout << p1->distance(*p2) << endl;
 auto start = system_clock::now();
 cout << p1->angle(*p2,*p3) << endl;
 auto end   = system_clock::now();
 auto duration = duration_cast<microseconds>(end - start);
 cout <<  "used "
 << double(duration.count()) * microseconds::period::num / microseconds::period::den
 << " second " << endl;
 */
