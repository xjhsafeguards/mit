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
    
    double OO_cutoff = 3.5;
    //double HCl_cutoff = 3;
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
            else if(strncmp(argv[i],"-st",3) == 0){string tmp=argv[++i];f_step=stoi(tmp);}
            else if(strncmp(argv[i],"-oo",3) == 0){string tmp=argv[++i];OO_cutoff=stod(tmp);}
            //else if(strncmp(argv[i],"-hcl",4) == 0){string tmp=argv[++i];HCl_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-h",2) == 0){cout << "-w for set water searching parameter OH_distance\n-p for set folder of position files\n-fs for the fisrt snapshot to read\n-fe for the last snapshot to read\n-st for the steps in read file\n-oo for OO_cutoff in reading OOO angle" << endl; return 1;}
            //else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];time_step=stod(tmp);}
            //else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    cout << "Reading folder: " + file_folder << endl;
    cout << "Reading snapshot from: " + to_string(f_start) + " to: " + to_string(f_end) + " with step: " + to_string(f_step) << endl;
    cout << "OO_cutoff :" + to_string(OO_cutoff) << endl;
    //cout << "HCl_cutoff :" + to_string(HCl_cutoff) << endl;
    
    double cell_vol=0;
    int scount=0;
    
    //auto Dp = new Distributionfunction(0,6,500);
    Distributionfunction* Dp[4];
    for(int i=0;i!=3;++i){
        Dp[i] = new Distributionfunction(0,6,500);
    }
    for(int i=3;i!=4;++i){
        Dp[i] = new Distributionfunction(0,180,500);
    }
    
    //cout << "test1" << endl;   
 
    for(int fc=0; fc!=8;++fc){
    
        ifstream ifs( file_folder + "/data.pos_"+ to_string(fc) + ".xyz");
        //molecule_manip* water = new water_manip();
        
        for(int i=0;i!=f_end;++i){
            
            std::shared_ptr<cell> cel = make_shared<cell_ipi>();
            cel->read(ifs);
            vector<vector<double>> data(4); // OO, OH, HH, OOO,
            //cout << "test1" << endl; 
            if(i>f_start and (i-f_start)%f_step == 0){
                
                for( const auto& atom1: cel->atoms()){
                    if(atom1->check_type("O"))
                        for( const auto& atom2: cel->atoms()){
                            if(atom2->check_type("O")){
                                data[0].push_back(atom1->distance(*atom2));
                                if(atom1!=atom2 and atom1->distance(*atom2)<OO_cutoff){
                                    //time consuming 
				    for( const auto& atom3: cel->atoms()){
                                        if(atom1!=atom3 and atom2!=atom3 and atom1->distance(*atom3)<OO_cutoff){
                                            data[3].push_back(atom1->angle(*atom2,*atom3));
                                   	    // cout << "test2" << endl;        
                                        }
                                    }
                                }
                            }
                            if(atom2->check_type("H"))
                                data[1].push_back(atom1->distance(*atom2));
                        }
                    if(atom1->check_type("H")){
                        for( const auto& atom2: cel->atoms()){
                            if(atom2->check_type("H")){
                                data[2].push_back(atom1->distance(*atom2));
                            }
                        }
                    }
                    cell_vol +=  cel->volume();
                    scount++;
                    for(int i=0;i!=4;++i){
                        //cout << data[i].size() << '\t';i
                        //cout << "test3" << endl;
                        Dp[i]->read(data[i]);
                        //cout << Dp[i]->get_valid_count() << '\t';
                    }
                }
                cout << "read bead" << fc << " snapshot " << i << '\r' << flush;
            }else{
                cout << "skip bead" << fc << " snapshot " << i << '\r' << flush;
            }
        }
    }
    
    cell_vol /= scount;
    
    vector<double> X=Dp[0]->get_x();
    vector<double> Y[4];
    for(int i=0;i!=3;++i){
        Dp[i]->set_dimension(3);
        Dp[i]->set_normalize(cell_vol);
        Y[i] = Dp[i]->get_y();
    }
    for(int i=3;i!=4;++i){
        //Dp[i]->set_dimension(3);
        //Dp[i]->set_normalize(cell_vol);
        Y[i] = Dp[i]->get_y();
    }
    
    ofstream ofs("G_OH_OO_OCl_HCl.txt");
    ofs << setprecision(10);
    
    ofs << "#" << setw(19) << "R" << setw(20) << "OO" << setw(20) << "OH" << setw(20) << "HH" << setw(20) << "OOO" <<  endl ;
    
    for(int i=0;i<X.size();++i){
        ofs << setw(20) << X[i];
        for(int j=0;j!=4;++j){
            ofs << setw(20) << Y[j][i];
        }
        ofs << '\n';
    }

    

}
