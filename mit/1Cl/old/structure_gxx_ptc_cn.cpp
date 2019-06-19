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
    
    Distributionfunction* Dp[6]; // G of OH, OO, OCl, HCl, OClO, HClH,
    for(int i=0;i!=4;++i){
        Dp[i] = new Distributionfunction(0,6,DF_step);
    }
    for(int i=4;i!=6;++i){
        Dp[i] = new Distributionfunction(0,180,DF_step);
    }
    Distributionfunction Dp_o(-4,1,DF_step); // PTC_O
    Distributionfunction Dp_cl(-4,1,DF_step); // PTC_Cl
    Distributionfunction Dp_CN(0.5,15.5,15); // CN
    
    molecule_manip* water = new water_manip();
    
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
                vector<vector<double>> data(6); // OH, OO, OCl, HCl, OClO, HClH,
                if(i>=f_start and (i-f_start)%f_step == 0){
                    water->read(*cel);
                    int CNcount=0;
                    // Reading GXX
                    GTIMER.start("GXX");
                    for( const auto& atom1: cel->atoms()){
                        if(atom1->check_type("O"))
                            for( const auto& atom2: cel->atoms()){
                                if(atom2->check_type("Cl"))
                                    data[2].push_back(atom1->distance(*atom2));
                                if(atom2->check_type("O")){
                                    data[1].push_back(atom1->distance(*atom2));
                                    if(atom1!=atom2 and atom1->distance(*cel->atoms()[0])<OCl_cutoff and atom2->distance(*cel->atoms()[0])<OCl_cutoff)
                                        data[4].push_back(cel->atoms()[0]->angle(*atom1,*atom2));
                                }
                                if(atom2->check_type("H"))
                                    data[0].push_back(atom1->distance(*atom2));
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
                    for(int i=0;i!=6;++i){
                        Dp[i]->read(data[i]);
                    }
                    GTIMER.stop("GXX");
                    // Reading PTC
                    GTIMER.start("PTC");
                    for(const auto& mol : cel->mols("H2O")){
                        //assert(mol->atoms().size()==3);
                        double OCl = mol->atoms()[0]->distance(*(cel->atoms()[0]));
                        double OH1 = mol->atoms()[1]->distance(*(mol->atoms()[0]));
                        double OH2 = mol->atoms()[2]->distance(*(mol->atoms()[0]));
                        if ( OCl < OCl_cutoff){
                            double HCl1 = mol->atoms()[1]->distance(*(cel->atoms()[0]));
                            double HCl2 = mol->atoms()[2]->distance(*(cel->atoms()[0]));
                            double PTC1 = OH1 - HCl1;
                            double PTC2 = OH2 - HCl2;
                            Dp_cl.read(PTC1);
                            Dp_cl.read(PTC2);
                            ++CNcount;
                        }
                        for(const auto& mol2 : cel->mols("H2O")){
                            double OO = mol->atoms()[0]->distance(*(mol2->atoms()[0]));
                            if( OO < OO_cutoff and OO > 0.001 ){
                                double H1O2= mol->atoms()[1]->distance(*(mol2->atoms()[0]));
                                double H2O2= mol->atoms()[2]->distance(*(mol2->atoms()[0]));
                                double PTC1 = OH1 - H1O2;
                                double PTC2 = OH2 - H2O2;
                                Dp_o.read(PTC1);
                                Dp_o.read(PTC2);
                            }
                            
                        }
                    }
                    GTIMER.stop("PTC");
                    Dp_CN.read(CNcount);
                    cell_vol +=  cel->volume();
                    scount++;
                    //cout << endl;
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
    
    
    cell_vol /= scount;
    { // print GXX
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
        ofstream ofs("G_OH_OO_OCl_HCl.txt");
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
    { // Print ptc
        ofstream ofs("ptc_cl_o.txt");
        ofs << setprecision(10);
        ofs << "#" << setw(19) << "delta" << setw(20) << "cl" << setw(20) << "o" << endl;
        ofs << "#OCl_cutoff: " << OCl_cutoff << endl;
        ofs << "#OO_cutoff: " << OO_cutoff << endl;
        
        auto X=Dp_cl.get_x();
        auto Y1=Dp_cl.get_y();
        auto Y2=Dp_o.get_y();
        
        for(int i=0; i<X.size(); ++i){
            ofs << setw(20) << X[i] << setw(20) << Y1[i] << setw(20) << Y2[i] << endl;
        }
        
    }
    { // Print CN
        vector<double> X=Dp_CN.get_x();
        vector<double> Y=Dp_CN.get_y();
        
        ofstream ofs("CN.txt");
        ofs << setprecision(10);
        
        ofs << "#" << setw(19) << "CN" << setw(20) << "Prob" << endl ;
        ofs << "#OCl_cutoff: " << OCl_cutoff << endl;
        
        for(int i=0;i<X.size();++i){
            ofs << setw(20) << X[i];
            ofs << setw(20) << Y[i];
            ofs << '\n';
        }
        
    }
    for(int i=0;i!=6;++i){
        Dp[i]->set_step_sigma(DF_sigma);
    }
    Dp_o.set_step_sigma(DF_sigma); // PTC_O
    Dp_cl.set_step_sigma(DF_sigma); // PTC_Cl
    { // print GXX
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
        ofstream ofs("G_OH_OO_OCl_HCl_sigma.txt");
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
    { // Print ptc
        ofstream ofs("ptc_cl_o_sigma.txt");
        ofs << setprecision(10);
        ofs << "#" << setw(19) << "delta" << setw(20) << "cl" << setw(20) << "o" << endl;
        ofs << "#OCl_cutoff: " << OCl_cutoff << endl;
        ofs << "#OO_cutoff: " << OO_cutoff << endl;
        
        auto X=Dp_cl.get_x();
        auto Y1=Dp_cl.get_y();
        auto Y2=Dp_o.get_y();
        
        for(int i=0; i<X.size(); ++i){
            ofs << setw(20) << X[i] << setw(20) << Y1[i] << setw(20) << Y2[i] << endl;
        }
        
    }
    GTIMER.stop("all");
    GTIMER.summerize(cout);
}

