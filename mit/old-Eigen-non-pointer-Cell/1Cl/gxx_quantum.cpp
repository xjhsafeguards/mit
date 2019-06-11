#include "Cell.h"
#include "Math_all.h"

#include <Utility_time.h>
#include <fstream>

using namespace std;

int main(int argc,char** argv){
    GTIMER.start("all");
    double OCl_cutoff = 3.5;
    double HCl_cutoff = 3;
    string file_folder = ".";
    int f_start = 0;
    int f_end = 10;
    
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            //if(strncmp(argv[i],"-n",2) == 0){string tmp=argv[++i];filename=tmp;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];snapshot_count=stoi(tmp);}
            //if(strncmp(argv[i],"-w",2) == 0){string tmp=argv[++i];water_parameter::OH_distance=stod(tmp);}
            if(strncmp(argv[i],"-p",2) == 0){file_folder=argv[++i];}
            else if(strncmp(argv[i],"-fs",3) == 0){string tmp=argv[++i];f_start=stoi(tmp);}
            else if(strncmp(argv[i],"-fe",3) == 0){string tmp=argv[++i];f_end=stoi(tmp);}
            else if(strncmp(argv[i],"-ocl",4) == 0){string tmp=argv[++i];OCl_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-hcl",4) == 0){string tmp=argv[++i];HCl_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-h",2) == 0){cout << "-w for set water searching parameter OH_distance\n-p for set folder of position files\n-fs for the fisrt snapshot to read\n-fe for the last snapshot to read" << endl; return 1;}
            //else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];time_step=stod(tmp);}
            //else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    cout << "Reading folder: " + file_folder << endl;
    cout << "Reading snapshot from: " + to_string(f_start) + " to: " + to_string(f_end) << endl;
    cout << "OCl_cutoff :" + to_string(OCl_cutoff) << endl;
    cout << "HCl_cutoff :" + to_string(HCl_cutoff) << endl;
    
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
        
        ifstream ifs( file_folder + "/data.pos_"+ to_string(fc) + ".xyz");
        //molecule_manip* water = new water_manip();
        
        for(int i=0;i!=f_end;++i){
            
            Cell_ipi cel;
            //std::shared_ptr<cell> cel = make_shared<cell_qecp>(cell_qecp({1,63,126},{"Cl","O","H"}));
            cel.read(ifs);
            vector<vector<double>> data(6); // OH, OO, OCl, HCl, OClO, HClH,
            
            if(i>f_start){
                
                for(int i=1; i<64; ++i){
                    //H
                    for(int j=64; j<190; ++j)
                        data[0].push_back(cel.adistance(i,j));
                    //Cl
                    double dOCl = cel.adistance(i,0);
                    data[2].push_back(dOCl);
                    //O
                    for(int j=1; j<64; ++j){
                        data[1].push_back(cel.adistance(i,j));
                        if(dOCl<OCl_cutoff and i!=j and cel.adistance(0,j)<OCl_cutoff)
                            data[4].push_back(cel.aangle(0,i,j));
                    }
                    
                }
                for(int i=64; i<190; ++i){
                    //Cl
                    double dHCl = cel.adistance(i,0);
                    data[3].push_back(dHCl);
                    //H
                    for(int j=64; j<190; ++j)
                        //data[1].push_back(cel.adistance(i,j));
                        if(dHCl<HCl_cutoff and i!=j and cel.adistance(0,j)<HCl_cutoff)
                            data[5].push_back(cel.aangle(0,i,j));
                }
                
                
                cell_vol +=  cel.volume();
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
    GTIMER.stop("all");
    GTIMER.summerize(cout);
}
