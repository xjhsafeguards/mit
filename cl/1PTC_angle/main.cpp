#include <iomanip>
#include <chrono>

#include "Math.h"
#include "Cell.h"
#include "Cell_Wannier90.h"
#include "Cell_TEMP.h"
#include "Molecule_water.h"

using namespace std;
using namespace chrono;


template<typename T>
void Print(const vector<T>& inv,ostream& os=cout){
    for( const auto& d: inv)
        os << setw(10) << d;
    os << endl;
}

//Read wannier file
void READ_WAN(){
    ifstream ifs("test_files/Ih_centres.xyz");
    
    cell* cel = new cell_wannier90();
    cel->set_box(12.444661365,12.444661365,12.444661365);
    cel->read(ifs);
    
    molecule_manip* mol = new water_manip();
    mol->read(*cel);
    mol->sort_wans(*cel);
    
    cout << cel->mols("H2O").size()*7 << endl << endl;
    for(const auto& mol : cel->mols("H2O")){
        for(const auto& atom : mol->atoms())
            cout << atom->get_type() << " " << atom->cart() << endl;
        for(const auto& wan : mol->wans())
            cout << wan->get_type() << " " << wan->cart() << endl;
        //cout << mol->atoms().size() << " " << mol->wans().size() << endl;
    }
}

int main(int argc,char** argv){
    
    string filename="test_files/cl.MD.snapshots_1219";
    int snapshot_count=21;
    
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            if(strncmp(argv[i],"-n",2) == 0){string tmp=argv[++i];filename=tmp;}
            else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];snapshot_count=stoi(tmp);}
            //else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];time_step=stod(tmp);}
            //else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    cout << setprecision(15);
    ifstream ifs(filename);
    
    cell* cel = new cell_temp1();
    cel->set_box(12.444661365,12.444661365,12.444661365);
    
    cout << setw(20) << "OH - ClH" << setw(20) << "O-H-Cl angle" << setw(20) << "O-Cl" << endl;
    
    for(int i=0;i!=snapshot_count;++i){
        cel->clear();
        //cout << "check 1" << endl;
        cel->read_atoms(ifs);
        //cout << "check 2" << endl;
        
        molecule_manip* mol = new water_manip();
        mol->read(*cel);
        //cout << "check 3" << endl;
        for(const auto& mol : cel->mols("H2O")){
            double HCl1 = mol->atoms()[1]->distance(*(cel->atoms()[0]));
            double HCl2 = mol->atoms()[2]->distance(*(cel->atoms()[0]));
            if( HCl1 < HCl2 )
                cout << setw(20) << mol->atoms()[1]->distance(*(mol->atoms()[0])) - HCl1 << setw(20) << mol->atoms()[1]->angle(*(mol->atoms()[0]),*(cel->atoms()[0])) << setw(20) << mol->atoms()[0]->distance(*(cel->atoms()[0])) << endl;
            else
                cout << setw(20) << mol->atoms()[2]->distance(*(mol->atoms()[0])) - HCl2 << setw(20) << mol->atoms()[2]->angle(*(mol->atoms()[0]),*(cel->atoms()[0])) << setw(20) << mol->atoms()[0]->distance(*(cel->atoms()[0])) << endl;
        }
    }

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
