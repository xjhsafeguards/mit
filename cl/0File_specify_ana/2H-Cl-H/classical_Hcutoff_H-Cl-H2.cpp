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
    
    Distributionfunction* Dp[15];
    for(int i=0;i!=15;++i){
        Dp[i] = new Distributionfunction(0,180,500);
        HCl_cutoff.push_back(2+0.2*i);
    }
    
    
    ifstream ifs1("/Users/jianhangxu/Documents/2Cl/Cl_63H2O_PBE_vdW_cpmd/cl.cel");
    ifstream ifs2("/Users/jianhangxu/Documents/2Cl/Cl_63H2O_PBE_vdW_cpmd/cl.pos");
        //molecule_manip* water = new water_manip();
        
    for(int i=0;i!=21220;++i){
        
        std::shared_ptr<cell> cel = make_shared<cell_qecp>(cell_qecp({1,63,126},{"Cl","O","H"}));
        cel->read_box(ifs1);
        cel->read_atoms(ifs2);
        if(i>5000){
            //water->read(*cel);
            //for(const auto& mol : cel->mols("H2O")){
            //assert(mol->atoms().size()==3);
            for( const auto& atom1: cel->atoms()){
                if(atom1->check_type("H") and atom1->distance(*cel->atoms()[0])<HCl_cutoff[14]){
                    for( const auto& atom2: cel->atoms()){
                        if(atom2->check_type("H") and atom2!=atom1 and atom2->distance(*cel->atoms()[0])<HCl_cutoff[14])
                        {
                            double DClH1 = atom1->distance(*cel->atoms()[0]);
                            double DClH2 = atom2->distance(*cel->atoms()[0]);
                            for(int j=0;j!=15;++j){
                                if(DClH1<HCl_cutoff[j] and DClH2<HCl_cutoff[j]){
                                    Dp[j]->read(cel->atoms()[0]->angle(*atom1,*atom2));
                                }
                            }
                        }
                    }
                    
                }
            }
        }
        cout << "Read snapshot " << i << '\r' << flush;
    }

    
    ofstream ofs("adf_Hcutoff.txt");
    ofs << setprecision(10);
    ofs << "#" << setw(19) << "angle";
    for(int i=0;i<15;++i){
        ofs << setw(20) << "Hcutoff:" + to_string(HCl_cutoff[i]);
    }
    ofs << endl << "#" << setw(19) << "Num count:";
    for(int i=0;i!=15;++i){
        ofs << setw(20) << Dp[i]->get_valid_count();
    }
    ofs << endl;
    
    vector<double> X=Dp[0]->get_x();
    vector<double> Y[15];
    for(int i=0;i!=15;++i){
        Y[i] = Dp[i]->get_y();
    }
    for(int i=0;i<X.size();++i){
        ofs << setw(20) << X[i];
        for(int j=0;j!=15;++j){
            ofs << setw(20) << Y[j][i];
        }
        ofs << '\n';
    }
}
    
