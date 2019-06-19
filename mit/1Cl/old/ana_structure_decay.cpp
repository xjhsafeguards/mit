#include <iomanip>
#include <chrono>
#include <cstring>

#include "Math_all.h"
#include "Cell.h"
#include "Cell_Wannier90.h"
#include "Cell_TEMP.h"
#include "Cell_QECP.h"
#include "Cell_IPI.h"
#include "Molecule_water.h"
#include "Utility_time.h"

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
    GTIMER.start("all");
    
    string Help_content;
    string file_folder = ".";
    int f_start = 0;
    int f_end = 10;
    int f_step = 1;
    double OCl_cutoff = 3.8;
    double HCl_cutoff = 2.9;
    double OO_cutoff = 3.35;
    int DF_step = 500;
    double DF_sigma = 0;
#ifdef QUANT
    water_parameter::OH_distance=2.0;
#endif
    Help_content += "-p for set folder of position files\n";
    Help_content += "-fs for the fisrt snapshot to read\n";
    Help_content += "-fe for the last snapshot to read\n";
    Help_content += "-st for the step of snapshot to read\n";
    Help_content += "-w for set water searching parameter OH_distance\n";
    Help_content += "-ocl for set the OCl_cutoff\n";
    Help_content += "-ohl for set the HCl_cutoff\n";
    Help_content += "-coo for set the OO_cutoff\n";
    Help_content += "-dfst for set the steps of distribution functions calculations\n";
    Help_content += "-dfsg for set the gausian sigma (in unit of step) of distribution functions calculations\n";
    
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            if(strncmp(argv[i],"-w",2) == 0){string tmp=argv[++i];water_parameter::OH_distance=stod(tmp);}
            else if(strncmp(argv[i],"-p",2) == 0){file_folder=argv[++i];}
            else if(strncmp(argv[i],"-fs",3) == 0){string tmp=argv[++i];f_start=stoi(tmp);}
            else if(strncmp(argv[i],"-fe",3) == 0){string tmp=argv[++i];f_end=stoi(tmp);}
            else if(strncmp(argv[i],"-st",3) == 0){string tmp=argv[++i];f_step=stoi(tmp);}
            else if(strncmp(argv[i],"-ocl",4) == 0){string tmp=argv[++i];OCl_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-hcl",4) == 0){string tmp=argv[++i];HCl_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-coo",4) == 0){string tmp=argv[++i];OO_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-dfst",5) == 0){string tmp=argv[++i];DF_step=stoi(tmp);}
            else if(strncmp(argv[i],"-dfsg",5) == 0){string tmp=argv[++i];DF_sigma=stod(tmp);}
            else if(strncmp(argv[i],"-h",2) == 0){cout << Help_content << endl; return 1;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    cout << "Reading folder: " + file_folder << endl;
    cout << "Reading snapshot from: " + to_string(f_start) + " to: " + to_string(f_end) + " with step: " + to_string(f_step) << endl;
    cout << "OCl_cutoff :" + to_string(OCl_cutoff) << endl;
    cout << "HCl_cutoff :" + to_string(HCl_cutoff) << endl;
    cout << "OO_cutoff :" + to_string(OO_cutoff) << endl;
    cout << "Distribution functions steps: " + to_string(DF_step) + " with gassian sigma(step): " + to_string(DF_sigma) << endl;
    
    double cell_vol=0;
    int scount=0;
    
    Distributionfunction* Dp[4];
    for(int i=0;i!=4;++i){
        Dp[i] = new Distributionfunction(0,180,DF_step); // OClO, ClOO, OOO, OO(cl-bonded)O
    }
    
#ifdef QUANT
    for(int fc=0; fc!=8;++fc){
        
        ifstream ifs( file_folder + "/data.pos_"+ to_string(fc) + ".xyz");
        //molecule_manip* water = new water_manip();
        
        for(int i=0;i!=f_end;++i){
            GTIMER.start("Read Cell");
            std::shared_ptr<cell> cel = make_shared<cell_ipi>();
            cel->read(ifs);
            GTIMER.stop("Read Cell");
#else
            ifstream ifs1(file_folder+"/cl.cel");
            ifstream ifs2(file_folder+"/cl.pos");
            
            for(int i=0;i!=f_end;++i){
                
                std::shared_ptr<cell> cel = make_shared<cell_qecp>(cell_qecp({1,63,126},{"Cl","O","H"}));
                GTIMER.start("Read Cell");
                cel->read_box(ifs1);
                cel->read_atoms(ifs2);
                GTIMER.stop("Read Cell");
#endif
                // Read Angle distribution
                vector<shared_ptr<position> > Os_in_Cl;
                for(int i=1;i<64;++i){
                    bool is_in=false;
                    if(cel->atoms()[0]->distance(*cel->atoms()[i])<OCl_cutoff){
                        Os_in_Cl.push_back(cel->atoms()[i]);
                        is_in=true;
                    }
                    vector<int> Oindex_in_Oi;
                    for(int j=1;j<64;++j){
                        if( i!=j and cel->atoms()[i]->distance(*cel->atoms()[j])<OO_cutoff)
                            Oindex_in_Oi.push_back(j);
                    }
                    for(int j: Oindex_in_Oi){
                        for(int k: Oindex_in_Oi)
                        {
                            if(j<k){
                                double Angle=cel->atoms()[i]->angle(*cel->atoms()[j],*cel->atoms()[k]);
                                Dp[2]->read(Angle); // OOO
                                if(is_in)
                                    Dp[3]->read(Angle);
                            }
                        }
                    }
                }
                for( const auto& atom1: Os_in_Cl){
                    for( const auto& atom2: Os_in_Cl){
                        if(atom1!=atom2 and atom1->distance(*atom2)<OO_cutoff)
                            Dp[0]->read(cel->atoms()[0]->angle(*atom1,*atom2)); // OClO
                    }
                    for(int i=1;i<64;++i){
                        if(atom1->distance(*cel->atoms()[i])<OO_cutoff)
                            Dp[1]->read(atom1->angle(*cel->atoms()[0],*cel->atoms()[i])); // ClOO
                    }
                }
                GTIMER.start("Print");
#ifdef QUANT
                cout << "read bead" << fc << " snapshot " << i << '\r' << flush;
                GTIMER.stop("Print");
            }
        }
#else
        
        cout << "read snapshot " << i << '\r' << flush;
        GTIMER.stop("Print");
    }
#endif
    
    vector<double> X=Dp[0]->get_x();
    vector<double> Y[4];
    for(int i=0;i!=4;++i){
        Y[i] = Dp[i]->get_y();
    }
    
    ofstream ofs("G_OClO_ClOO_OOO.txt");
    ofs << setprecision(10);
    
    ofs << "#" << setw(19) << "R" << setw(20) << "OClO" << setw(20) << "ClOO" << setw(20) << "OOO" << setw(20) << "OO(cl)O" << endl;
    
    for(int i=0;i<X.size();++i){
        ofs << setw(20) << X[i];
        for(int j=0;j!=4;++j){
            ofs << setw(20) << Y[j][i];
        }
        ofs << '\n';
    }
    GTIMER.stop("all");
    GTIMER.summerize(cout);
}

