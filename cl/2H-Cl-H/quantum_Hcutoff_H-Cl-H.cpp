#include <iomanip>
#include <chrono>

#include "Math.h"
#include "Cell.h"
#include "Cell_Wannier90.h"
#include "Cell_TEMP.h"
#include "Cell_QECP.h"
#include "Cell_IPI.h"
#include "Molecule_water.h"

using namespace std;
using namespace chrono;


template<typename T>
void Print(const vector<T>& inv,ostream& os=cout){
    for( const auto& d: inv)
        //os << setw(10) << d;
        os << d << endl;
    os << endl;
}


int main(int argc,char** argv){
    
    double OCl_cutoff = 4;
    vector<double> HCl_cutoff;
    

    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            //if(strncmp(argv[i],"-n",2) == 0){string tmp=argv[++i];filename=tmp;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];snapshot_count=stoi(tmp);}
            if(strncmp(argv[i],"-w",2) == 0){string tmp=argv[++i];water_parameter::OH_distance=stod(tmp);}
            //else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];time_step=stod(tmp);}
            //else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    Distributionfunction* Dp[7];
    for(int i=0;i!=7;++i){
        Dp[i] = new Distributionfunction(0,180,500);
        HCl_cutoff.push_back(2.7+0.1*i);
    }
    
    for(int fc=0; fc!=8;++fc){
    
        ifstream ifs("/Users/jianhangxu/Documents/2Cl/Cl_63H2O_pimd/data.pos_"+ to_string(fc) + ".xyz");
        //molecule_manip* water = new water_manip();
        
        for(int i=0;i!=48230;++i){
            
            std::shared_ptr<cell> cel = make_shared<cell_ipi>();
            cel->read(ifs);
            if(i>10000){
                //water->read(*cel);
                //for(const auto& mol : cel->mols("H2O")){
                //assert(mol->atoms().size()==3);
                for( const auto& atom1: cel->atoms()){
                    if(atom1->check_type("H") and atom1->distance(*cel->atoms()[0])<HCl_cutoff[6]){
                        for( const auto& atom2: cel->atoms()){
                            if(atom2->check_type("H") and atom2!=atom1 and atom2->distance(*cel->atoms()[0])<HCl_cutoff[6])
                            {
                                double DClH1 = atom1->distance(*cel->atoms()[0]);
                                double DClH2 = atom2->distance(*cel->atoms()[0]);
                                    for(int j=0;j!=7;++j){
                                        if(DClH1<HCl_cutoff[j] and DClH2<HCl_cutoff[j]){
                                            Dp[j]->read(cel->atoms()[0]->angle(*atom1,*atom2));
                                        }
                                    }
                            }
                        }

                    }
                }
            }
            cout << "read bead" << fc << " snapshot " << i << '\r' << flush;
        }
    }
    
    ofstream ofs("adf_Hcutoff.txt");
    ofs << setprecision(10);
    ofs << "#" << setw(19) << "angle";
    for(int i=0;i<7;++i){
        ofs << setw(20) << "Hcutoff:" + to_string(HCl_cutoff[i]);
    }
    ofs << endl << "#" << setw(19) << "Num count:";
    for(int i=0;i!=7;++i){
        ofs << setw(20) << Dp[i]->get_valid_count();
    }
    ofs << endl;
    
    vector<double> X=Dp[0]->get_x();
    vector<double> Y[10];
    for(int i=0;i!=7;++i){
        Y[i] = Dp[i]->get_y();
    }
    for(int i=0;i<X.size();++i){
        ofs << setw(20) << X[i];
        for(int j=0;j!=7;++j){
            ofs << setw(20) << Y[j][i];
        }
        ofs << '\n';
    }
    
    /*
    double cell_vol=0;
    int scount=0;
    vector<double> OH;
    vector<double> OO;
    vector<double> OCl;
    vector<double> OClO;
    vector<double> HClH;
    vector<double> HCl;
    for(auto it = cellv.begin();it!=cellv.end();++it){
        for( const auto& atom1: (*it)->atoms()){
            if(atom1->check_type("O"))
                for( const auto& atom2: (*it)->atoms()){
                    if(atom2->check_type("H"))
                        OH.push_back(atom1->distance(*atom2));
                    if(atom2->check_type("Cl"))
                        OCl.push_back(atom1->distance(*atom2));
                    if(atom2->check_type("O")){
                        OO.push_back(atom1->distance(*atom2));
                        if(atom1->distance(*(*it)->atoms()[0])<OCl_cutoff and atom2->distance(*(*it)->atoms()[0])<OCl_cutoff)
                            OClO.push_back((*it)->atoms()[0]->angle(*atom1,*atom2));
                    }
                }
            if(atom1->check_type("H")){
                for( const auto& atom2: (*it)->atoms()){
                    if(atom2->check_type("H") and atom1->distance(*(*it)->atoms()[0])<HCl_cutoff and atom2->distance(*(*it)->atoms()[0])<HCl_cutoff)
                        HClH.push_back((*it)->atoms()[0]->angle(*atom1,*atom2));
                    if(atom2->check_type("Cl"))
                        HCl.push_back(atom1->distance(*atom2));
                }
            }
        }
        cout << "read rdf " << ++scount << '\r' << flush;
        cell_vol +=  (*it)->volume();
    }
    
    cell_vol /= scount;
    
    auto Dp = new Distributionfunction(0,6,500);
    Dp->set_dimension(3);
    Dp->set_normalize(cell_vol);
    Dp->read(OH);
    
    vector<double> X=Dp->get_x();
    vector<double> YOH=Dp->get_y();

    Dp->reset();
    Dp->read(OO);
    vector<double> YOO=Dp->get_y();
    
    Dp->reset();
    Dp->read(OCl);
    vector<double> YOCl=Dp->get_y();
    
    Dp->reset();
    Dp->read(HCl);
    vector<double> YHCl=Dp->get_y();
    
    ofstream ofs("rdf.txt");
    ofs << setprecision(10);
    ofs << setw(20) << "#X" << setw(20) << "OH" << setw(20) << "OO" << setw(20) << "OCl" << setw(20) << "HCl" << endl;
    for(int i=0;i<X.size();++i){
        ofs << setw(20) << X[i] << setw(20) << YOH[i] << setw(20) << YOO[i] << setw(20) << YOCl[i] << setw(20) << YHCl[i] << endl;
    }
    
    ofstream ofs2("adf.txt");
    ofs << setprecision(10);
    
    Dp->reset(0,180,400);
    Dp->set_dimension(1);
    Dp->set_normalize(1);
    Dp->read(OClO);
    
    X=Dp->get_x();
    vector<double> YOClO=Dp->get_y();
    
    Dp->reset();
    Dp->read(HClH);
    vector<double> YHClH=Dp->get_y();
    
    ofs2 << setw(20) << "#X" << setw(20) << "OClO" << setw(20) << "HClH" << endl;
    for(int i=0;i<X.size();++i){
        ofs2 << setw(20) << X[i] << setw(20) << YOClO[i] << setw(20) << YHClH[i] << endl;
    }
    */
}

//Read wannier file
void READ_WAN1(){
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

void READ_WAN2(int argc,char** argv){
    string filename="test_files/cl.MD.snapshots_1219";
    int snapshot_count=21;
    
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            if(strncmp(argv[i],"-n",2) == 0){string tmp=argv[++i];filename=tmp;}
            else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];snapshot_count=stoi(tmp);}
            else if(strncmp(argv[i],"-w",2) == 0){string tmp=argv[++i];water_parameter::OH_distance=stod(tmp);}
            //else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];time_step=stod(tmp);}
            //else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    cout << setprecision(15);
    ifstream ifs(filename);
    
    cout << setw(20) << "#O-H" << setw(20) << "O-W" << setw(20) << "Cl-W" << setw(20) << "O-Cl" << setw(10) << ""<< endl;
    
    cell* cel = new cell_wannier90();
    cel->set_box(12.444661365,12.444661365,12.444661365);
    
    for(int i=0;i!=snapshot_count;++i){
        cel->clear();
        cel->read(ifs);
        
        molecule_manip* mol = new water_manip();
        mol->read(*cel);
        mol->sort_wans(*cel);
        
        cout << endl;
        //cout << cel->mols("H2O").size()*7 << endl << endl;
        for(const auto& mol : cel->mols("H2O")){
            
            double HCl1 = mol->atoms()[1]->distance(*(cel->atoms()[0]));
            double HCl2 = mol->atoms()[2]->distance(*(cel->atoms()[0]));
            if( HCl1 < HCl2 )
                cout << setw(20) << mol->atoms()[1]->distance(*(mol->atoms()[0])) << setw(20) << mol->atoms()[0]->distance(*(mol->wans()[0])) << setw(20) << mol->wans()[0]->distance(*(cel->atoms()[0])) << setw(20) << mol->atoms()[0]->distance(*(cel->atoms()[0])) << setw(10) << mol->atoms().size()-1 << endl;
            else
                cout << setw(20) << mol->atoms()[2]->distance(*(mol->atoms()[0])) << setw(20) << mol->atoms()[0]->distance(*(mol->wans()[1])) << setw(20) << mol->wans()[1]->distance(*(cel->atoms()[0])) << setw(20) << mol->atoms()[0]->distance(*(cel->atoms()[0])) << setw(10) << mol->atoms().size()-1 << endl;
            //for(const auto& atom : mol->atoms())
            //cout << atom->get_type() << " " << atom->cart() << endl;
            //for(const auto& wan : mol->wans())
            //cout << wan->get_type() << " " << wan->cart() << endl;
            //cout << mol->atoms().size() << " " << mol->wans().size() << endl;
        }
    }
}


void READ_PTC_ANGLE(int argc,char** argv){
    string filename="test_files/cl.MD.snapshots_1219";
    int snapshot_count=21;
    
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            if(strncmp(argv[i],"-n",2) == 0){string tmp=argv[++i];filename=tmp;}
            else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];snapshot_count=stoi(tmp);}
            else if(strncmp(argv[i],"-w",2) == 0){string tmp=argv[++i];water_parameter::OH_distance=stod(tmp);}
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
    
    cout << setw(20) << "#OH - ClH" << setw(20) << "O-H-Cl angle" << setw(20) << "O-Cl" << endl;
    
    for(int i=0;i!=snapshot_count;++i){
        cel->clear();
        //cout << "check 1" << endl;
        cel->read_atoms(ifs);
        //cout << "check 2" << endl;
        
        molecule_manip* mol = new water_manip();
        mol->read(*cel);
        //cout << "check 3" << endl;
        cout << endl;
        for(const auto& mol : cel->mols("H2O")){
            //assert(mol->atoms().size()==3);
            double HCl1 = mol->atoms()[1]->distance(*(cel->atoms()[0]));
            double HCl2 = mol->atoms()[2]->distance(*(cel->atoms()[0]));
            if( HCl1 < HCl2 )
                cout << setw(20) << mol->atoms()[1]->distance(*(mol->atoms()[0])) - HCl1 << setw(20) << mol->atoms()[1]->angle(*(mol->atoms()[0]),*(cel->atoms()[0])) << setw(20) << mol->atoms()[0]->distance(*(cel->atoms()[0])) << setw(10) << mol->atoms().size()-1 << endl;
            else
                cout << setw(20) << mol->atoms()[2]->distance(*(mol->atoms()[0])) - HCl2 << setw(20) << mol->atoms()[2]->angle(*(mol->atoms()[0]),*(cel->atoms()[0])) << setw(20) << mol->atoms()[0]->distance(*(cel->atoms()[0])) << setw(10) << mol->atoms().size()-1 << endl;
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
