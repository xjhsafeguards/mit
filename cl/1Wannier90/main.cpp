#include <iomanip>
#include <chrono>

#include "Math.h"
#include "Cell.h"
#include "Cell_Wannier90.h"
#include "Cell_TEMP.h"
#include "Molecule_water.h"

using namespace std;
int main(int argc,char** argv){
    
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

