#include <iomanip>
#include <chrono>
#include <cstring>

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
    
    double OCl_cutoff = 3.5;
    double HCl_cutoff = 3;
    double PTC_cutoff = -0.5;
    string file_folder = ".";
    int f_start = 0;
    int f_end = 10;
    
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
            else if(strncmp(argv[i],"-ocl",4) == 0){string tmp=argv[++i];OCl_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-hcl",4) == 0){string tmp=argv[++i];HCl_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-cptc",5) == 0){string tmp=argv[++i];PTC_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-h",2) == 0){cout << "-w for set water searching parameter OH_distance\n-ocl for set the OCl_cutoff\n-ohl for set the HCl_cutoff\n-p for set folder of position files\n-fs for the fisrt snapshot to read\n-fe for the last snapshot to read\n-cptc for PTC_cutoff" << endl; return 1;}
            //else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];time_step=stod(tmp);}
            //else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    cout << "Reading folder: " + file_folder << endl;
    cout << "Reading snapshot from: " + to_string(f_start) + " to: " + to_string(f_end) << endl;
    cout << "OCl_cutoff :" + to_string(OCl_cutoff) << endl;
    cout << "HCl_cutoff :" + to_string(HCl_cutoff) << endl;
    
    //auto Dp = new Distributionfunction(-4,0,500);
    molecule_manip* water = new water_manip();
    int CN,nptc=20,ncn=15;
    double PTCAVG,PTC_low=-OCl_cutoff,PTC_up=0;
    vector<vector<int> > PTC_CN(nptc,vector<int>(ncn));

    //Distributionfunction Dp(-4,0,500);
    
    for(int fc=0; fc!=8;++fc){
        
        ifstream ifs( file_folder + "/data.pos_"+ to_string(fc) + ".xyz");
        
        for(int i=0;i!=f_end;++i){
            
            std::shared_ptr<cell> cel = make_shared<cell_ipi>();
            cel->read(ifs);
            if(i>f_start){
                water->read(*cel);
                CN=0;
                PTCAVG=0;
                for(const auto& mol : cel->mols("H2O")){
                    //assert(mol->atoms().size()==3);
                    double OCl = mol->atoms()[0]->distance(*(cel->atoms()[0]));
                    if ( OCl < OCl_cutoff){
                        double HCl1 = mol->atoms()[1]->distance(*(cel->atoms()[0]));
                        double HCl2 = mol->atoms()[2]->distance(*(cel->atoms()[0]));
                        double PTC = (HCl1 < HCl2) ? mol->atoms()[1]->distance(*(mol->atoms()[0])) - HCl1 : mol->atoms()[2]->distance(*(mol->atoms()[0])) - HCl2;
                        PTCAVG+=PTC;
                        ++CN;
                    }
                }
                ++PTC_CN[static_cast<int>((PTCAVG/CN-PTC_low)/(PTC_up-PTC_low)*nptc)][CN];
            }
            cout << "read bead" << fc << " snapshot " << i << '\r' << flush;
        }
        //     cellv.push_back(cel);
        //cout << "read " << cellv.size() << '\r' << flush;
    }
    
    ofstream ofs("ptc_cn.txt");
    ofs << setprecision(10);
    ofs << "#" << setw(19) << "PTC" << setw(20) << "CN" << setw(20) << "Count" << endl;
    for(int i=0; i<PTC_CN.size(); ++i){
        for(int j=0; j<PTC_CN[i].size(); ++j)
        {
            ofs << setw(20) << static_cast<double>(i)/nptc*(PTC_up-PTC_low)+PTC_low << setw(20) << j << setw(20) << PTC_CN[i][j] << endl;
        }
    }
}

    
