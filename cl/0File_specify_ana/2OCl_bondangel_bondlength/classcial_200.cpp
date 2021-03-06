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
    
    //vector<double> OCl_cutoff;
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
    
    //auto Dp = new Distributionfunction(-4,0,500);
    molecule_manip* water = new water_manip();
    //Distributionfunction* Dp[6];
    //for(int i=0;i!=6;++i){
    //    Dp[i] = new Distributionfunction(-4,0,500);
    //    OCl_cutoff.push_back(3.5+0.1*i);
    //}
    Stat_pool ba[200];
    Stat_pool bl[200];
    
    ifstream ifs1("/Users/jianhangxu/Documents/2Cl/Cl_63H2O_PBE_vdW_cpmd/cl.cel");
    ifstream ifs2("/Users/jianhangxu/Documents/2Cl/Cl_63H2O_PBE_vdW_cpmd/cl.pos");
    
    for(int i=0;i!=21220;++i){
        
        std::shared_ptr<cell> cel = make_shared<cell_qecp>(cell_qecp({1,63,126},{"Cl","O","H"}));
        cel->read_box(ifs1);
        cel->read_atoms(ifs2);
        cout << "Read snapshot " << i << '\r' << flush;
        if(i>5000){
            water->read(*cel);
            for(const auto& mol : cel->mols("H2O")){
                //assert(mol->atoms().size()==3);
                double OCl = mol->atoms()[0]->distance(*(cel->atoms()[0]));
                int O_range = (OCl-2)/0.05;
                if( O_range>-1 and O_range<200){
                    ba[O_range].read(mol->atoms()[0]->angle(*mol->atoms()[1],*mol->atoms()[2]));
                    bl[O_range].read((mol->atoms()[0]->distance(*mol->atoms()[1])+mol->atoms()[0]->distance(*mol->atoms()[2]))/2);                    }
            }
        }
    }
    
    ofstream ofs("H2Odistance_Bondangle_Bondlength.txt");
    ofs << setprecision(10);
    ofs << "#" << setw(19) << "O-Cl" << setw(20) << "Num count:" << setw(20) << "bond_angle_ave" << setw(20) << "bond_angle_sd" << setw(20) << "bond_angle_max" << setw(20) << "bond_angle_min" << setw(20) << "bond_length_ave" << setw(20) << "bond_length_sd" << setw(20) << "bond_length_max" << setw(20) << "bond_length_min\n";
    

    for(int i=0;i!=200;++i){
        if( ba[i].count() < 1000)
            continue;
        ofs << setw(20) << 2.025+0.05*i;
        ofs << setw(20) << ba[i].count();
        ofs << setw(20) << ba[i].ave();
        ofs << setw(20) << ba[i].sd();
        ofs << setw(20) << ba[i].max();
        ofs << setw(20) << ba[i].min();
        ofs << setw(20) << bl[i].ave();
        ofs << setw(20) << bl[i].sd();
        ofs << setw(20) << bl[i].max();
        ofs << setw(20) << bl[i].min();
        ofs << '\n';
    }
}
    
