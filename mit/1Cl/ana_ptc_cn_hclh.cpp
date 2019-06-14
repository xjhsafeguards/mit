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
    double PTC_cutoff_up = -1;
    double PTC_cutoff_down = -1.2;
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
    Help_content += "-dptcu for set upper limit of PTC\n";
    Help_content += "-dptcd for set lower limit of PTC\n";
    
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
            else if(strncmp(argv[i],"-dptcu",6) == 0){string tmp=argv[++i];PTC_cutoff_up=stod(tmp);}
            else if(strncmp(argv[i],"-dptcd",6) == 0){string tmp=argv[++i];PTC_cutoff_down=stod(tmp);}
            else if(strncmp(argv[i],"-h",2) == 0){cout << Help_content << endl; return 1;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    cout << "Reading folder: " + file_folder << endl;
    cout << "Reading snapshot from: " + to_string(f_start) + " to: " + to_string(f_end) + " with step: " + to_string(f_step) << endl;
    cout << "OCl_cutoff :" + to_string(OCl_cutoff) << endl;
    cout << "HCl_cutoff :" + to_string(HCl_cutoff) << endl;
    cout << "OO_cutoff :" + to_string(OO_cutoff) << endl;
    cout << "Distribution functions steps: " + to_string(DF_step) << endl; //+ " with gassian sigma(step): " + to_string(DF_sigma) << endl;
    
    Distributionfunction* Dp_nH[15]; //HClH of CNH
    Distributionfunction* Dp_nO[15]; //HClH of CNO
    Distributionfunction* Dp_nPTC[15]; //HClH of snapshot with n of required PTC
    for(int i=0;i!=15;++i){
        Dp_nH[i] = new Distributionfunction(0,180,DF_step);
    }
    for(int i=0;i!=15;++i){
        Dp_nO[i] = new Distributionfunction(0,180,DF_step);
    }
    for(int i=0;i!=15;++i){
        Dp_nPTC[i] = new Distributionfunction(0,180,DF_step);
    }
    Distributionfunction Dp_one(0,180,DF_step); //HClH of atoms only with required PTC
    Distributionfunction Dp_both(0,180,DF_step); //HClH of atoms both have required PTC
    Distributionfunction Dp(0,180,DF_step);
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
                GTIMER.start("Read Water");
                water->read(*cel);
                GTIMER.stop("Read Water");
                
                int CNcount = 0,CNHcount = 0,nPTC = 0;
                vector<double>  HClH;
                vector<shared_ptr<position> > Hs_in_Cl,Hs_with_ptc;
                
                GTIMER.start("Read HClH");
                for(int i=64;i<190;++i){
                    bool is_in=false;
                    double dHCl = cel->atoms()[0]->distance(*cel->atoms()[i]);
                    if( dHCl<HCl_cutoff){
                        Hs_in_Cl.push_back(cel->atoms()[i]);
                        ++CNHcount;
                    }
                }
                for(const auto& mol : cel->mols("H2O")){
                    //assert(mol->atoms().size()==3);
                    double OCl = mol->atoms()[0]->distance(*(cel->atoms()[0]));
                    if ( OCl < OCl_cutoff){
                        ++CNcount;
                        double OH1 = mol->atoms()[1]->distance(*(mol->atoms()[0]));
                        double OH2 = mol->atoms()[2]->distance(*(mol->atoms()[0]));
                        double HCl1 = mol->atoms()[1]->distance(*(cel->atoms()[0]));
                        double HCl2 = mol->atoms()[2]->distance(*(cel->atoms()[0]));
                        double PTC1 = OH1 - HCl1;
                        double PTC2 = OH2 - HCl2;
                        if(PTC1 < PTC_cutoff_up and PTC1 > PTC_cutoff_down){
                            Hs_with_ptc.push_back(mol->atoms()[1]);
                            ++nPTC;
                        }
                        if(PTC2 < PTC_cutoff_up and PTC2 > PTC_cutoff_down){
                            Hs_with_ptc.push_back(mol->atoms()[2]);
                            ++nPTC;
                        }
                    }
                }
                for( const auto& atom1: Hs_in_Cl){
                    for( const auto& atom2: Hs_in_Cl)
                        if(atom1!=atom2)
                            HClH.push_back(cel->atoms()[0]->angle(*atom1,*atom2)); // HClH
                    for( const auto& atom2: Hs_with_ptc)
                        if(atom1!=atom2)
                            Dp_one.read(cel->atoms()[0]->angle(*atom1,*atom2)); // HClH
                }
                for( const auto& atom1: Hs_with_ptc)
                    for( const auto& atom2: Hs_with_ptc)
                        if(atom1!=atom2)
                            Dp_both.read(cel->atoms()[0]->angle(*atom1,*atom2)); // HClH
                
                Dp_nH[CNHcount]->read(HClH);
                Dp_nPTC[nPTC]->read(HClH);
                Dp_nO[CNcount]->read(HClH);
                Dp.read(HClH);
                GTIMER.stop("Read HClH");
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
    
    {// print adf HClH CN
        ofstream ofs("adf_nH.txt");
        ofs << setprecision(10);
        ofs << "#" << setw(19) << "angle" << setw(20) << "total";
        for(int i=0;i<15;++i){
            ofs << setw(20) << "nH:" + to_string(i);
        }
        for(int i=0;i<15;++i){
            ofs << setw(20) << "nO:" + to_string(i);
        }
        ofs << endl << "#" << setw(19) << "Num count:" << setw(20) << Dp.get_valid_count() ;
        for(int i=0;i!=15;++i){
            ofs << setw(20) << Dp_nH[i]->get_valid_count();
        }
        for(int i=0;i!=15;++i){
            ofs << setw(20) << Dp_nO[i]->get_valid_count();
        }
        ofs << endl;
        
        vector<double> X=Dp_nH[0]->get_x();
        vector<double> Y=Dp.get_y();
        vector<double> Y1[15],Y2[15];
        for(int i=0;i!=15;++i){
            Y1[i] = Dp_nH[i]->get_y();
            Y2[i] = Dp_nO[i]->get_y();
        }
        for(int i=0;i<X.size();++i){
            ofs << setw(20) << X[i];
            ofs << setw(20) << Y[i];
            for(int j=0;j!=15;++j){
                ofs << setw(20) << Y1[j][i];
            }
            for(int j=0;j!=15;++j){
                ofs << setw(20) << Y2[j][i];
            }
            ofs << '\n';
        }
    }
    {// print adf with in PTC
        string ofname;
        ofname = "PTC" + to_string(PTC_cutoff_up) + " to " + to_string(PTC_cutoff_down) + "H-Cl-H_one.txt";
        ofstream ofs(ofname);
        ofs << setprecision(10);
        ofs << "#" << setw(19) << "angle" << setw(20) << "one" << setw(20) << "both";
        for(int i=0;i<15;++i){
            ofs << setw(20) << "nvalidptc:" + to_string(i);
        }
        ofs << endl << "#" << setw(19) << "Num count:";
        ofs << setw(20) << Dp_one.get_valid_count() << setw(20) << Dp_both.get_valid_count();
        for(int i=0;i!=15;++i){
            ofs << setw(20) << Dp_nPTC[i]->get_valid_count();
        }
        ofs << endl;
        
        vector<double> X=Dp.get_x();
        vector<double> Y1=Dp_one.get_y();
        vector<double> Y2=Dp_both.get_y();
        vector<double> Y3[15];
        for(int i=0;i!=15;++i){
            Y3[i] = Dp_nPTC[i]->get_y();
        }
        for(int i=0;i<X.size();++i){
            ofs << setw(20) << X[i];
            ofs << setw(20) << Y1[i];
            ofs << setw(20) << Y2[i];
            for(int j=0;j!=15;++j){
                ofs << setw(20) << Y3[j][i];
            }
            ofs << '\n';
        }
    }
    GTIMER.stop("all");
    GTIMER.summerize(cout);
}




