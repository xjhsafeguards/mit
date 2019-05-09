#include <iomanip>
#include <chrono>
#include <cstring>
#include <algorithm>
#include <utility>

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
        os << d << "\n";
    os << endl;
}


int main(int argc,char** argv){
    
    double OCl_cutoff = 3.8;
    double HCl_cutoff = 2.9;
    double OO_cutoff = 3.35;
    //double PTC_cutoff = -0.5;
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
            else if(strncmp(argv[i],"-h",2) == 0){cout << "-w for set water searching parameter OH_distance\n-ocl for set the OCl_cutoff\n-ohl for set the HCl_cutoff\n-coo for set the OO_cutoff\n-p for set folder of position files\n-fs for the fisrt snapshot to read\n-fe for the last snapshot to read\n-ft for the steps in reading\n" << endl; return 1;}
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
    
    //auto Dp = new Distributionfunction2D(-OCl_cutoff,0,500,-0.5,14.5,15);
    molecule_manip* water = new water_manip();
    //int CN,nptc=20,ncn=15;
    //double PTCAVG,PTC_low=-OCl_cutoff,PTC_up=0;
    //vector<vector<int> > PTC_CN(nptc,vector<int>(ncn));
    
    //Distributionfunction Dp_o(-4,1,500);
    //Distributionfunction Dp_cl(-4,1,500);
    
    ofstream ofs[1];
    for(int i=0;i<1;++i){
        ofs[i].open("CN_6_"+to_string(i)+".xyz");
        ofs[i] << setprecision(10);
    }
    vector<Vector3<double> > ref_pos; //= {Vector3<double>(HCl_cutoff,0,0),Vector3<double>(0,HCl_cutoff,0),Vector3<double>(-HCl_cutoff,0,0),Vector3<double>(0,-HCl_cutoff,0),Vector3<double>(0,0,HCl_cutoff),Vector3<double>(0,0,-HCl_cutoff)};
    
    for(int fc=0; fc!=8;++fc){
        
        ifstream ifs( file_folder + "/data.pos_"+ to_string(fc) + ".xyz");
        for(int i=0;i!=f_end;++i){
            
            std::shared_ptr<cell> cel = make_shared<cell_ipi>();
            cel->read(ifs);
            
            if(i>f_start and (i-f_start)%f_step == 0){
                
                water->read(*cel);
                vector<pair<double,shared_ptr<position> > > certified_H;
                vector<pair<double,shared_ptr<position> > > certified_O;
                for(const auto& mol : cel->mols("H2O")){
                    double HCl1 = mol->atoms()[1]->distance(*(cel->atoms()[0]));
                    double HCl2 = mol->atoms()[2]->distance(*(cel->atoms()[0]));
                    double OCl = mol->atoms()[0]->distance(*(cel->atoms()[0]));
                    double OH1 = mol->atoms()[1]->distance(*(mol->atoms()[0]));
                    double OH2 = mol->atoms()[2]->distance(*(mol->atoms()[0]));
                    double PTC1 = OH1 - HCl1;
                    double PTC2 = OH2 - HCl2;
                    if( HCl1 < HCl_cutoff){
                        certified_H.push_back(make_pair(PTC1,mol->atoms()[1]));
                        certified_O.push_back(make_pair(PTC1,mol->atoms()[0]));
                    }
                    if( HCl2 < HCl_cutoff){
                        certified_H.push_back(make_pair(PTC2,mol->atoms()[2]));
                        certified_O.push_back(make_pair(PTC2,mol->atoms()[0]));
                    }
                }
                /*
                 for(const auto& atom : cel->atoms()){
                 if(atom->check_type("H") and atom->distance(*cel->atoms()[0]) < HCl_cutoff ){
                 certified_H.push_back(atom);
                 }
                 }*/
                int count = certified_H.size();
                if(count == 6){
                    //sort by PTC
                    sort(certified_H.begin(),certified_H.end(),[](pair<double,shared_ptr<position> >& a,pair<double,shared_ptr<position> >& b){return a.first>b.first;});
                    sort(certified_O.begin(),certified_O.end(),[](pair<double,shared_ptr<position> >& a,pair<double,shared_ptr<position> >& b){return a.first>b.first;});
                    // construct a list of 6 positions
                    vector<Vector3<double> > pos_list_H,pos_list_O;
                    for(const auto& atom : certified_H){
                        pos_list_H.push_back(atom.second->cart().shortest_BC(cel->atoms()[0]->cart(),cel->boxp()->diagonal()));
                    }
                    for(const auto& atom : certified_O){
                        pos_list_O.push_back(atom.second->cart().shortest_BC(cel->atoms()[0]->cart(),cel->boxp()->diagonal()));
                    }
                    //sort(pos_list.begin(),pos_list.end(),[](Vector3<double>& a,Vector3<double>& b){return a.norm()<b.norm();});
                    // rotate the position of PTC min to 001
                    Optimal_rotation orot1;
                    orot1.Solve(vector<Vector3<double>>{pos_list_O[0]},vector<Vector3<double>>{Vector3<double>(0,0,1)});
                    auto rot_pos_H = orot1.Rotated_vectors(pos_list_H);
                    auto rot_pos_O = orot1.Rotated_vectors(pos_list_O);
                    
                    // if never store a reference create one
                    if( ref_pos.size() == 0 )
                        ref_pos = rot_pos_O;
                    //permutation
                    //construct a 6*6 matrix
                    vector<vector<double> > permut_m(6,vector<double>(6));
                    for(int i=0; i<6; ++i)
                        for(int j=0; j<6; ++j)
                            permut_m[i][j] = pow(ref_pos[i].distance(rot_pos_O[j]),2);
                    //solve the linear assignment problem
                    vector<int> per_sol;
                    HungarianAlgorithm ha;
                    ha.Solve(permut_m,per_sol);
                    vector<Vector3<double> > sorted_pos_O(6),sorted_pos_H(6);
                    for(int i=0;i<6;++i){
                        //sorted_pos[per_sol[i]] = rot_pos[i];
                        sorted_pos_O[i] = rot_pos_O[per_sol[i]];
                        sorted_pos_H[i] = rot_pos_H[per_sol[i]];
                    }
                    
                    //solve the rotation problem
                    Optimal_rotation orot;
                    orot.allow_reflection();
                    auto result_pos_O = orot.Solve(sorted_pos_O,ref_pos);
                    auto result_pos_H = orot.Rotated_vectors(sorted_pos_H);
                    
                    //auto rotation_m = orot.Rotation_matrix();
                    /*
                     cout << endl;
                     double d1=0,d2=0,d3=0;
                     for(int i=0;i<6;++i){
                     //cout << result_pos[i] << " | " << sorted_pos[i] << " | " << ref_pos[i] << endl;
                     d1 +=  pow(result_pos[i].distance(ref_pos[i]),2);
                     d2 +=  pow(sorted_pos[i].distance(ref_pos[i]),2);
                     d3 +=  pow(ref_pos[i].distance(ref_pos[i]),2);
                     }
                     cout << d1 <<  " | " << d2 << " | " << d3;
                     //if(d1 > d2) cout << "   ******************";
                     cout << endl;
                     */
                    ofs[0] << count*2+1 << '\n' << "SS: " << i << " PTC: ";
                    for(const auto& item: certified_H)
                        ofs[0] << item.first << " ";
                    ofs[0] <<'\n';
                    ofs[0] << " Cl 0 0 0" << '\n';
                    for(const auto& atom : result_pos_H){
                        ofs[0] << " H  ";
                        ofs[0] << atom << '\n';
                    }
                    for(const auto& atom : result_pos_O){
                        ofs[0] << " O  ";
                        ofs[0] << atom << '\n';
                    }
                }
            }
            cout << "read bead" << fc << " snapshot " << i << '\r' << flush;
        }
    }
}



