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
    
    double cell_vol=0;
    int scount=0;
    
    //auto Dp = new Distributionfunction(0,6,500);
    Distributionfunction* Dp[6];
    for(int i=0;i!=4;++i){
        Dp[i] = new Distributionfunction(0,6,500);
    }
    for(int i=4;i!=6;++i){
        Dp[i] = new Distributionfunction(0,180,500);
    }
    
    
    for(int fc=0; fc!=8;++fc){
    
        ifstream ifs("/Users/jianhangxu/Documents/2Cl/cl_63H2O_npt_pi_300K/9.2ps/data.pos_"+ to_string(fc) + ".xyz");
        //molecule_manip* water = new water_manip();
        
        for(int i=0;i!=19300;++i){
            
            std::shared_ptr<cell> cel = make_shared<cell_ipi>();
            cel->read(ifs);
            vector<vector<double>> data(6); // OH, OO, OCl, HCl, OClO, HClH,
            
            if(i>12000){
                
                for( const auto& atom1: cel->atoms()){
                    if(atom1->check_type("O"))
                        for( const auto& atom2: cel->atoms()){
                            if(atom2->check_type("H"))
                                data[0].push_back(atom1->distance(*atom2));
                            if(atom2->check_type("Cl"))
                                data[2].push_back(atom1->distance(*atom2));
                            if(atom2->check_type("O")){
                                data[1].push_back(atom1->distance(*atom2));
                                if(atom1!=atom2 and atom1->distance(*cel->atoms()[0])<OCl_cutoff and atom2->distance(*cel->atoms()[0])<OCl_cutoff)
                                    data[4].push_back(cel->atoms()[0]->angle(*atom1,*atom2));
                            }
                        }
                    if(atom1->check_type("H")){
                        for( const auto& atom2: cel->atoms()){
                            if(atom2->check_type("H") and atom1!=atom2 and atom1->distance(*cel->atoms()[0])<HCl_cutoff and atom2->distance(*cel->atoms()[0])<HCl_cutoff)
                                data[5].push_back(cel->atoms()[0]->angle(*atom1,*atom2));
                            if(atom2->check_type("Cl"))
                                data[3].push_back(atom1->distance(*atom2));
                        }
                    }
                }
                cell_vol +=  cel->volume();
                scount++;
                for(int i=0;i!=6;++i){
                    //cout << data[i].size() << '\t';
                    Dp[i]->read(data[i]);
                    //cout << Dp[i]->get_valid_count() << '\t';
                }
                //cout << endl;
            }
            cout << "read bead" << fc << " snapshot " << i << '\r' << flush;
        }
    }
    
    cell_vol /= scount;
    
    vector<double> X=Dp[0]->get_x();
    vector<double> Y[6];
    for(int i=0;i!=4;++i){
        Dp[i]->set_dimension(3);
        Dp[i]->set_normalize(cell_vol);
        Y[i] = Dp[i]->get_y();
    }
    for(int i=4;i!=6;++i){
        //Dp[i]->set_dimension(3);
        //Dp[i]->set_normalize(cell_vol);
        Y[i] = Dp[i]->get_y();
    }
    
    ofstream ofs("G_OH_OO_OCl_HCl_OClO_HClH.txt");
    ofs << setprecision(10);
    
    ofs << "#" << setw(19) << "R" << setw(20) << "OH" << setw(20) << "OO" << setw(20) << "OCl" << setw(20) << "HCl" << setw(20) << "OClO" << setw(20) << "HClH" << endl ;
    
    for(int i=0;i<X.size();++i){
        ofs << setw(20) << X[i];
        for(int j=0;j!=6;++j){
            ofs << setw(20) << Y[j][i];
        }
        ofs << '\n';
    }

    

}
