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
    
    Distributionfunction* Dp[9]; // G of OH, OO, OCl, HCl, OClO, HClH, ClOO, OOO, OO(cl-bonded)O
    for(int i=0;i!=4;++i){
        Dp[i] = new Distributionfunction(0,6,DF_step);
    }
    for(int i=4;i!=9;++i){
        Dp[i] = new Distributionfunction(0,180,DF_step);
    }
    Distributionfunction Dp_o(-4,1,DF_step); // PTC_O
    Distributionfunction Dp_cl(-4,1,DF_step); // PTC_Cl
    Distributionfunction Dp_CN(0.5,15.5,15); // CN
    Distributionfunction Dp_CNH(0.5,15.5,15); // CN
    Distributionfunction Dp_CNd(-7.5,7.5,15); // CN diff O-H
    vector<Stat_pool> ba(DF_step),bl(DF_step);
    double decay_step = (6.0-0)/DF_step;
    
    molecule_manip* water = new water_manip();
    
#ifdef QUANT
    for(int fc=0; fc!=8;++fc){
        
        ifstream ifs( file_folder + "/data.pos_"+ to_string(fc) + ".xyz");
        //molecule_manip* water = new water_manip();
        
        for(int i=0;i!=f_end;++i){
            std::shared_ptr<cell> cel = make_shared<cell_ipi>();
            if(i>=f_start and (i-f_start)%f_step == 0){
                GTIMER.start("Read Cell");
                cel->read(ifs);
                GTIMER.stop("Read Cell");
#else
                ifstream ifs1(file_folder+"/cl.cel");
                ifstream ifs2(file_folder+"/cl.pos");
                
                for(int i=0;i!=f_end;++i){
                    std::shared_ptr<cell> cel = make_shared<cell_qecp>(cell_qecp({1,63,126},{"Cl","O","H"}));
                    if(i>=f_start and (i-f_start)%f_step == 0){
                        GTIMER.start("Read Cell");
                        cel->read_box(ifs1);
                        cel->read_atoms(ifs2);
                        GTIMER.stop("Read Cell");
#endif
                        GTIMER.start("Read Water");
                        water->read(*cel);
                        GTIMER.stop("Read Water");
                        int CNcount=0,CNHcount=0;
                        // Reading GXX assuming all atoms are in 1cl 63O 126H
                        GTIMER.start("GXX");
                        vector<shared_ptr<position> > Os_in_Cl;
                        vector<shared_ptr<position> > Hs_in_Cl;
                        for(int i=1;i<64;++i){
                            bool is_in=false;
                            double dOCl = cel->atoms()[0]->distance(*cel->atoms()[i]);
                            Dp[2]->read(dOCl); // OCl
                            if( dOCl<OCl_cutoff){
                                Os_in_Cl.push_back(cel->atoms()[i]);
                                is_in=true;
                            }
                            vector<int> Oindex_in_Oi;
                            for(int j=1;j<64;++j){
                                double dOO = cel->atoms()[i]->distance(*cel->atoms()[j]);
                                Dp[1]->read(dOO); // OO
                                if( i!=j and dOO < OO_cutoff)
                                    Oindex_in_Oi.push_back(j);
                            }
                            for(int j: Oindex_in_Oi){
                                for(int k: Oindex_in_Oi)
                                {
                                    if(j<k){
                                        double Angle=cel->atoms()[i]->angle(*cel->atoms()[j],*cel->atoms()[k]);
                                        Dp[7]->read(Angle); // OOO
                                        if(is_in)
                                            Dp[8]->read(Angle); // OO(cl)O
                                    }
                                }
                            }
                            for(int j=64;j<190;++j)
                                Dp[0]->read(cel->atoms()[i]->distance(*cel->atoms()[j])); // OH
                        }
                        for( const auto& atom1: Os_in_Cl){
                            for( const auto& atom2: Os_in_Cl){
                                if(atom1!=atom2)// and atom1->distance(*atom2)<OO_cutoff
                                    Dp[4]->read(cel->atoms()[0]->angle(*atom1,*atom2)); // OClO
                            }
                            for(int i=1;i<64;++i){
                                if(atom1!=cel->atoms()[i] and atom1->distance(*cel->atoms()[i])<OO_cutoff)
                                    Dp[6]->read(atom1->angle(*cel->atoms()[0],*cel->atoms()[i])); // ClOO
                            }
                        }
                        for(int i=64;i<190;++i){
                            bool is_in=false;
                            double dHCl = cel->atoms()[0]->distance(*cel->atoms()[i]);
                            Dp[3]->read(dHCl); // HCl
                            if( dHCl<HCl_cutoff){
                                Hs_in_Cl.push_back(cel->atoms()[i]);
                                ++CNHcount;
                            }
                            
                        }
                        for( const auto& atom1: Hs_in_Cl)
                            for( const auto& atom2: Hs_in_Cl)
                                if(atom1!=atom2)// and atom1->distance(*atom2)<OO_cutoff
                                    Dp[5]->read(cel->atoms()[0]->angle(*atom1,*atom2)); // HClH
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
                            int O_range = OCl/decay_step;
                            if( O_range>-1 and O_range<DF_step){
                                ba[O_range].read(mol->atoms()[0]->angle(*mol->atoms()[1],*mol->atoms()[2]));
                                bl[O_range].read((mol->atoms()[0]->distance(*mol->atoms()[1])+mol->atoms()[0]->distance(*mol->atoms()[2]))/2);                    }
                        }
                        GTIMER.stop("PTC");
                        Dp_CN.read(CNcount);
                        Dp_CNH.read(CNHcount);
                        Dp_CNd.read(CNcount-CNHcount);
                        cell_vol +=  cel->volume();
                        scount++;
                        //cout << endl;
                        GTIMER.start("Print");
#ifdef QUANT
                        cout << "read bead" << fc << " snapshot " << i << '\r' << flush;
                        GTIMER.stop("Print");
                    }
                    else{
                        GTIMER.start("Read Cell");
                        cel->skip(ifs);
                        GTIMER.stop("Read Cell");
                        GTIMER.start("Print");
                        cout << "skip bead" << fc << " snapshot " << i << '\r' << flush;
                        GTIMER.stop("Print");
                    }
                }
            }
#else
            
            cout << "read snapshot " << i << '\r' << flush;
            GTIMER.stop("Print");
        }
        else{
            GTIMER.start("Read Cell");
            cel->skip_box(ifs1);
            cel->skip_atoms(ifs2);
            GTIMER.stop("Read Cell");
            GTIMER.start("Print");
            cout << "skip snapshot " << i << '\r' << flush;
            GTIMER.stop("Print");
        }
    }
#endif
    
    
    cell_vol /= scount;
    { // print GXX
        vector<double> X=Dp[0]->get_x();
        vector<double> Y[9];
        for(int i=0;i!=4;++i){
            Dp[i]->set_dimension(3);
            Dp[i]->set_normalize(cell_vol);
            Y[i] = Dp[i]->get_y();
        }
        for(int i=4;i!=9;++i){
            //Dp[i]->set_dimension(3);
            //Dp[i]->set_normalize(cell_vol);
            Y[i] = Dp[i]->get_y();
        }
        ofstream ofs("G_OH_OO_OCl_HCl.txt");
        ofs << setprecision(10);
        
        ofs << "#" << setw(19) << "R" << setw(20) << "OH" << setw(20) << "OO" << setw(20) << "OCl" << setw(20) << "HCl" << setw(20) << "OClO" << setw(20) << "HClH" << setw(20) << "ClOO" << setw(20) << "OOO" << setw(20) << "OO(cl-bonded)O" << endl ;
        
        for(int i=0;i<X.size();++i){
            ofs << setw(20) << X[i];
            for(int j=0;j!=9;++j){
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
        vector<double> X2=Dp_CNd.get_x();
        vector<double> Y=Dp_CN.get_y();
        vector<double> Y2=Dp_CNH.get_y();
        vector<double> Y3=Dp_CNd.get_y();
        
        ofstream ofs("CN.txt");
        ofs << setprecision(10);
        
        ofs << "#" << setw(19) << "CN" << setw(20) << "OProb" << setw(20) << "HProb" << setw(20) << "CN diff(O-H)" << setw(20) << "Prob" <<endl ;
        ofs << "#OCl_cutoff: " << OCl_cutoff << "   HCl_cutoff: " << HCl_cutoff << endl;
        
        for(int i=0;i<X.size();++i){
            ofs << setw(20) << X[i];
            ofs << setw(20) << Y[i];
            ofs << setw(20) << Y2[i];
            ofs << setw(20) << X2[i];
            ofs << setw(20) << Y3[i];
            ofs << '\n';
        }
        
    }
    for(int i=0;i!=9;++i){
        Dp[i]->set_step_sigma(DF_sigma);
    }
    Dp_o.set_step_sigma(DF_sigma); // PTC_O
    Dp_cl.set_step_sigma(DF_sigma); // PTC_Cl
    { // print GXX
        vector<double> X=Dp[0]->get_x();
        vector<double> Y[9];
        for(int i=0;i!=4;++i){
            Dp[i]->set_dimension(3);
            Dp[i]->set_normalize(cell_vol);
            Y[i] = Dp[i]->get_y();
        }
        for(int i=4;i!=9;++i){
            //Dp[i]->set_dimension(3);
            //Dp[i]->set_normalize(cell_vol);
            Y[i] = Dp[i]->get_y();
        }
        ofstream ofs("G_OH_OO_OCl_HCl_sigma.txt");
        ofs << setprecision(10);
        
        ofs << "#" << setw(19) << "R" << setw(20) << "OH" << setw(20) << "OO" << setw(20) << "OCl" << setw(20) << "HCl" << setw(20) << "OClO" << setw(20) << "HClH" << setw(20) << "ClOO" << setw(20) << "OOO" << setw(20) << "OO(cl-bonded)O" << endl ;
        
        for(int i=0;i<X.size();++i){
            ofs << setw(20) << X[i];
            for(int j=0;j!=9;++j){
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
    { // print bond-length dacay
        ofstream ofs("H2Odistance_Bondangle_Bondlength.txt");
        ofs << setprecision(10);
        ofs << "#" << setw(19) << "O-Cl" << setw(20) << "Num count:" << setw(20) << "bond_angle_ave" << setw(20) << "bond_angle_sd" << setw(20) << "bond_angle_max" << setw(20) << "bond_angle_min" << setw(20) << "bond_length_ave" << setw(20) << "bond_length_sd" << setw(20) << "bond_length_max" << setw(20) << "bond_length_min\n";
        
        
        for(int i=0;i!=DF_step;++i){
            //if( ba[i].count() < 1000)
            //continue;
            ofs << setw(20) << decay_step*(i+0.5);
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
    GTIMER.stop("all");
    GTIMER.summerize(cout);
}

/*
 save for later
 for( const auto& atom: cel->atom_type("O")){
 double dOCl=atom->distance(*Cl);
 data[2].push_back(dOCl)); // OCl
 if(dOCl<OCl_cutoff)
 Os_in_Cl.push_back(atom);
 }
 for( const auto& atom: cel->atom_type("H")){
 double dHCl=atom->distance(*Cl);
 data[3].push_back(dHCl)); // HCl
 if(dHCl<OCl_cutoff)
 Hs_in_Cl.push_back(atom);
 }
 for( const auto& atom: cel->atom_type("O")){
 */




