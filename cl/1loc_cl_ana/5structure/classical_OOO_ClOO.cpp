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
    
    double OCl_cutoff = 3.8;
    double HCl_cutoff = 2.9;
    double OO_cutoff = 3.35;
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
            else if(strncmp(argv[i],"-h",2) == 0){cout << "-w for set water searching parameter OH_distance\n-ocl for set the OCl_cutoff\n-ohl for set the HCl_cutoff\n-coo for set the OO_cutoff\n-p for set folder of position files\n-fs for the fisrt snapshot to read\n-fe for the last snapshot to read\n-ft for the steps in reading" << endl; return 1;}
            //else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];time_step=stod(tmp);}
            //else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    cout << "Reading folder: " + file_folder << endl;
    cout << "Reading snapshot from: " + to_string(f_start) + " to: " + to_string(f_end) + " with step: " + to_string(f_step) << endl;
    cout << "OCl_cutoff :" + to_string(OCl_cutoff) << endl;
    cout << "HCl_cutoff :" + to_string(HCl_cutoff) << endl;
    cout << "OO_cutoff :" + to_string(OO_cutoff) << endl;
    
    //double cell_vol=0;
    //int scount=0;
    
    //auto Dp = new Distributionfunction(0,6,500);
    Distributionfunction* Dp[3];
    for(int i=0;i!=3;++i){
        //Dp[i] = new Distributionfunction(0,6,500);
        Dp[i] = new Distributionfunction(0,180,500); // OClO, ClOO, OOO
    }
    //for(int i=4;i!=6;++i){
       // Dp[i] = new Distributionfunction(0,180,500);
    //}
    
    
    ifstream ifs1(file_folder+"/cl.cel");
    ifstream ifs2(file_folder+"/cl.pos");
        
    for(int i=0;i!=f_end;++i){
        
        std::shared_ptr<cell> cel = make_shared<cell_qecp>(cell_qecp({1,63,126},{"Cl","O","H"}));
        cel->read_box(ifs1);
        cel->read_atoms(ifs2);
        
        if(i>f_start and (i-f_start)%f_step == 0){
            
            vector<shared_ptr<position> > Os_in_Cl;
            for(int i=1;i<64;++i){
                if(cel->atoms()[0]->distance(*cel->atoms()[i])<OCl_cutoff)
                    Os_in_Cl.push_back(cel->atoms()[i]);
                vector<int> Oindex_in_Oi;
                for(int j=1;j<64;++j){
                    if( i!=j and cel->atoms()[i]->distance(*cel->atoms()[j])<OO_cutoff)
                        Oindex_in_Oi.push_back(j);
                }
                for(int j: Oindex_in_Oi){
                    for(int k: Oindex_in_Oi)
                    {
                        if(j<k)
                            Dp[2]->read(cel->atoms()[i]->angle(*cel->atoms()[j],*cel->atoms()[k])); // OOO
                    }
                }
            }
            for( const auto& atom1: Os_in_Cl){
                for( const auto& atom2: Os_in_Cl){
                    if(atom1!=atom2 and atom1->distance(*atom2)<OO_cutoff)
                        //data[0].push_back(cel->atoms()[0]->angle(*atom1,*atom2)); // OClO
                        Dp[0]->read(cel->atoms()[0]->angle(*atom1,*atom2)); // OClO
                }
                for(int i=1;i<64;++i){
                    if(atom1->distance(*cel->atoms()[i])<OO_cutoff)
                        //data[1].push_back(atom1->angle(*cel->atoms()[0],*cel->atoms()[i])); // ClOO
                        Dp[1]->read(atom1->angle(*cel->atoms()[0],*cel->atoms()[i])); // ClOO
                }
            }
        }
        cout << "Read snapshot " << i << '\r' << flush;
    }
    
    //cell_vol /= scount;
    
    vector<double> X=Dp[0]->get_x();
    vector<double> Y[3];
    for(int i=0;i!=3;++i){
        //Dp[i]->set_dimension(3);
        //Dp[i]->set_normalize(cell_vol);
        Y[i] = Dp[i]->get_y();
    }
    //for(int i=4;i!=6;++i){
        //Dp[i]->set_dimension(3);
        //Dp[i]->set_normalize(cell_vol);
        //Y[i] = Dp[i]->get_y();
    //}
    
    ofstream ofs("G_OClO_ClOO_OOO.txt");
    ofs << setprecision(10);
    
    ofs << "#" << setw(19) << "R" << setw(20) << "OClO" << setw(20) << "ClOO" << setw(20) << "OOO" << endl ;
    
    for(int i=0;i<X.size();++i){
        ofs << setw(20) << X[i];
        for(int j=0;j!=3;++j){
            ofs << setw(20) << Y[j][i];
        }
        ofs << '\n';
    }

    

}
