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
    
    double OCl_cutoff = 3.82;
    double HCl_cutoff = 3;
    double PTC_cutoff = -1.2;
    
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
    
    auto Dp = new Distributionfunction(0,180,501);
    auto Dp2 = new Distributionfunction(0,180,501);
    molecule_manip* water = new water_manip();
    //Distributionfunction* Dp[6];
    //for(int i=0;i!=6;++i){
    //    Dp[i] = new Distributionfunction(-4,0,500);
    //    OCl_cutoff.push_back(3.5+0.1*i);
    //}
    
    ifstream ifs1("/Users/jianhangxu/Documents/2Cl/cl_63H2O_npt_bo_300K/out.15.134/cl.cel");
    ifstream ifs2("/Users/jianhangxu/Documents/2Cl/cl_63H2O_npt_bo_300K/out.15.134/cl.pos");
    //molecule_manip* water = new water_manip();
    
    for(int i=0;i!=31200;++i){
        
        std::shared_ptr<cell> cel = make_shared<cell_qecp>(cell_qecp({1,63,126},{"Cl","O","H"}));
        cel->read_box(ifs1);
        cel->read_atoms(ifs2);
        
        if(i>12000){
            water->read(*cel);
            int j=0;
            for(const auto& mol : cel->mols("H2O")){
                //assert(mol->atoms().size()==3);
                double OCl = mol->atoms()[0]->distance(*(cel->atoms()[0]));
                if ( OCl < OCl_cutoff){
                    double HCl1 = mol->atoms()[1]->distance(*(cel->atoms()[0]));
                    double HCl2 = mol->atoms()[2]->distance(*(cel->atoms()[0]));
                    double PTC = (HCl1 < HCl2) ? mol->atoms()[1]->distance(*(mol->atoms()[0])) - HCl1 : mol->atoms()[2]->distance(*(mol->atoms()[0])) - HCl2;
                    if(PTC<PTC_cutoff){
                        int Hlabel = (HCl1 < HCl2) ? 1 : 2;
                        for(const auto& mol2 : cel->mols("H2O")){
                            //assert(mol->atoms().size()==3);
                            if(mol!=mol2){
                                double OCl = mol2->atoms()[0]->distance(*(cel->atoms()[0]));
                                if ( OCl < OCl_cutoff){
                                    double HCl1 = mol2->atoms()[1]->distance(*(cel->atoms()[0]));
                                    double HCl2 = mol2->atoms()[2]->distance(*(cel->atoms()[0]));
                                    //double PTC = (HCl1 < HCl2) ? mol2->atoms()[1]->distance(*(mol2->atoms()[0])) - HCl1 : mol2->atoms()[2]->distance(*(mol2->atoms()[0])) - HCl2;
                                    if(PTC<PTC_cutoff){
                                        int Hlabel2 = (HCl1 < HCl2) ? 1 : 2;
                                        Dp->read(cel->atoms()[0]->angle(*(mol->atoms()[Hlabel]),*(mol2->atoms()[Hlabel2])));
                                    }
                                }
                            }
                        }
                        /*
                         for( const auto& atom2: cel->atoms()){
                         if(atom2->check_type("H") and atom2!=mol->atoms()[Hlabel] and atom2->distance(*cel->atoms()[0])<HCl_cutoff)
                         Dp->read(cel->atoms()[0]->angle(*(mol->atoms()[Hlabel]),*atom2));
                         }*/
                    }
                }
            }
        }
        //     cellv.push_back(cel);
        //cout << "read " << cellv.size() << '\r' << flush;
        cout << "Read snapshot " << i << '\r' << flush;
    }
    
    ofstream ofs("PTC>>" + to_string(PTC_cutoff) + "H-Cl-H.txt");
    ofs << setprecision(10);
    ofs << "#" << setw(19) << "angle" << setw(20) << "distribution";
    ofs << endl << "#" << setw(19) << "Num count:";
    ofs << setw(20) << Dp->get_valid_count();
    ofs << endl;
    
    vector<double> X=Dp->get_x();
    vector<double> Y=Dp->get_y();
    for(int i=0;i<X.size();++i){
        ofs << setw(20) << X[i];
        ofs << setw(20) << Y[i];
        ofs << '\n';
    }
}
    
