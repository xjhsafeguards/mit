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
    double HCl_cutoff = 3;

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
    
    vector<double> cellsize[8];
    
    
    for(int fc=0; fc!=8;++fc){
    
        ifstream ifs("/Users/jianhangxu/Documents/2Cl/cl_63H2O_npt_pi_300K/pi_bp-2019-03-11/data.pos_"+ to_string(fc) + ".xyz");
        //molecule_manip* water = new water_manip();
        
        for(int i=0;i!=27100;++i){
            
            std::shared_ptr<cell> cel = make_shared<cell_ipi>();
            cel->read(ifs);

            cellsize[fc].push_back(cel->volume());
            cout << "read bead" << fc << " snapshot " << i << '\r' << flush;
        }
    }
    
    ofstream ofs("SS_volume.txt");
    ofs << setprecision(10);
    
    ofs << "#" << setw(19) << "SS";
    for(int fc=0; fc!=8;++fc){
        ofs << setw(20) << "Bead" + to_string(fc);
    }
    ofs << endl ;
    
    for(int i=0;i<cellsize[0].size();++i){
        ofs << setw(20) << i;
        for(int j=0;j!=8;++j){
            ofs << setw(20) << cellsize[j][i];
        }
        ofs << '\n';
    }

    

}
