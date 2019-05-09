#include <iomanip>
#include <chrono>
#include <cstring>
#include <algorithm>
#include <utility>

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
        os << d << "\n";
    os << endl;
}


int main(int argc,char** argv){
    
    double OCl_cutoff = 3.8;
    double HCl_cutoff = 2.9;
    double OO_cutoff = 3.35;
    //double PTC_cutoff = -0.5;
    string file_folder = ".";
    int f_start = 0;
    int f_end = 10;
    int f_step = 1;
    
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            //if(strncmp(argv[i],"-n",2) == 0){string tmp=argv[++i];filename=tmp;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];snapshot_count=stoi(tmp);}
            if(strncmp(argv[i],"-w",2) == 0){string tmp=argv[++i];water_parameter::OH_distance=stod(tmp);}
            else if(strncmp(argv[i],"-p",2) == 0){file_folder=argv[++i];}
            else if(strncmp(argv[i],"-fs",3) == 0){string tmp=argv[++i];f_start=stoi(tmp);}
            else if(strncmp(argv[i],"-fe",3) == 0){string tmp=argv[++i];f_end=stoi(tmp);}
            else if(strncmp(argv[i],"-ft",3) == 0){string tmp=argv[++i];f_step=stoi(tmp);}
            else if(strncmp(argv[i],"-ocl",4) == 0){string tmp=argv[++i];OCl_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-hcl",4) == 0){string tmp=argv[++i];HCl_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-coo",4) == 0){string tmp=argv[++i];OO_cutoff=stod(tmp);}
            //else if(strncmp(argv[i],"-cptc",5) == 0){string tmp=argv[++i];PTC_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-h",2) == 0){cout << "-w for set water searching parameter OH_distance\n-ocl for set the OCl_cutoff\n-ohl for set the HCl_cutoff\n-coo for set the OO_cutoff\n-p for set folder of position files\n-fs for the fisrt snapshot to read\n-fe for the last snapshot to read\n-ft for the steps in reading\n" << endl; return 1;}
            //else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];time_step=stod(tmp);}
            //else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    cout << "Reading folder: " + file_folder << endl;
    cout << "Reading snapshot from: " + to_string(f_start) + " to: " + to_string(f_end) + " with step: " + to_string(f_step) << endl;
    //cout << "OCl_cutoff :" + to_string(OCl_cutoff) << endl;
    //cout << "HCl_cutoff :" + to_string(HCl_cutoff) << endl;
    //cout << "OO_cutoff :" + to_string(OO_cutoff) << endl;
    
    //auto Dp = new Distributionfunction2D(-OCl_cutoff,0,500,-0.5,14.5,15);
    //molecule_manip* water = new water_manip();
    //int CN,nptc=20,ncn=15;
    //double PTCAVG,PTC_low=-OCl_cutoff,PTC_up=0;
    //vector<vector<int> > PTC_CN(nptc,vector<int>(ncn));
    
    //Distributionfunction Dp_o(-4,1,500);
    //Distributionfunction Dp_cl(-4,1,500);
    
    ofstream ofs[1];
    for(int i=0;i<1;++i){
        ofs[i].open("Cl_SCAN_MD.snapshot");
        ofs[i] << setprecision(10);
    }
    
    ifstream ifs1(file_folder+"/cl.cel");
    ifstream ifs2(file_folder+"/cl.pos");
    
    for(int i=0;i!=f_end;++i){
        
        std::shared_ptr<cell> cel = make_shared<cell_qecp>(cell_qecp({1,63,126},{"Cl","O","H"}));
        cel->read_box(ifs1);
        cel->read_atoms(ifs2);
        
        if(i>f_start and (i-f_start)%f_step == 0){
            
            vector<pair<double,shared_ptr<position> > > O;
            vector<shared_ptr<position> > H;
            
            for(const auto& atom : cel->atoms()){
                if(atom->check_type("O")){
                    O.push_back(make_pair(atom->distance(*(cel->atoms()[0])),atom));
                }else if(atom->check_type("H")){
                    H.push_back(atom);
                }
            }
            sort(O.begin(),O.end(),[](pair<double,shared_ptr<position> >& a,pair<double,shared_ptr<position> >& b){return a.first<b.first;});
            
            ofs[0] << "timestep: " << i << '\n';
            ofs[0] << "Cl " << cel->atoms()[0]->normal_cart() << '\n';
            for(const auto& item : O){
                ofs[0] << "O  " << item.second->normal_cart() << "  " << item.first << '\n';
            }
            for(const auto& atom : H){
                ofs[0] << "H  " << atom->normal_cart() << '\n';
            }
            
        }
        cout << "Read snapshot " << i << '\r' << flush;
    }
    
}



