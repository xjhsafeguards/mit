#include "Cell.h"
#include "Cell_ipi.h"
#include "Cell_qe.h"
#include "../cl/Math.h"
#include "Molecule_water.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

using namespace std;

int main(int argc, char** argv){
    double OCl_cutoff = 3.5;
    double HCl_cutoff = 3;
    double PTC_cutoff_up = -1;
    double PTC_cutoff_down = -1.2;
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
            else if(strncmp(argv[i],"-ocl",4) == 0){string tmp=argv[++i];OCl_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-hcl",4) == 0){string tmp=argv[++i];HCl_cutoff=stod(tmp);}
            else if(strncmp(argv[i],"-pu",3) == 0){string tmp=argv[++i];PTC_cutoff_up=stod(tmp);}
            else if(strncmp(argv[i],"-pd",3) == 0){string tmp=argv[++i];PTC_cutoff_down=stod(tmp);}
            else if(strncmp(argv[i],"-h",2) == 0){cout << "-w for set water searching parameter OH_distance\n-p for set folder of position files\n-fs for the fisrt snapshot to read\n-fe for the last snapshot to read-fs for the fisrt snapshot to read\n-pu PTC_cutoff upper limit\n -pd PTC_cutoff lower limit" << endl; return 1;}
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
    cout << "PTC_cutoff_up :" + to_string(PTC_cutoff_up) << endl;
    cout << "PTC_cutoff_down :" + to_string(PTC_cutoff_down) << endl;
    
    auto Dp = new Distributionfunction(0,180,501);

    for(int fc=0; fc!=8;++fc){
        
        ifstream ifs( file_folder + "/data.pos_"+ to_string(fc) + ".xyz");
        bool have_seed=false;
        typename Water_group::seed_type WS;
        //molecule_manip* water = new water_manip();
    
        for(int i=0;i!=f_end;++i){
            
            Cell_ipi_c cel;
            
            if(i>f_start and (i-f_start)%f_step == 0){
                cel.read(ifs);
                cout << "test7" << endl ;
                cout << "test8" << endl ;
                Water_group wa(cel);
                if(!have_seed){
                    wa.read();
                    have_seed=true;
                    WS = wa.seed();
                }
                else{
                    wa=std::move(Water_group(cel,WS)); 
                }
                auto Cl=cel.atom(0);
                cout << "test2" << endl ;
                for(auto const& wp: wa){
                    cout << "test6 indexO water: " << wp->atom(0)->get_index() << endl ;
                    double OCl=wp->atom(0)->distance(Cl);
                    cout << "test3" << endl ;
                    if ( OCl < OCl_cutoff){
                        cout << "test7" << endl ;
                        double HCl1 = wp->atom(1)->distance(Cl);
                        double HCl2 = wp->atom(2)->distance(Cl);
                        cout << "test5" << endl ;
                        double PTC = (HCl1 < HCl2) ? wp->atom(1)->distance(*(wp->atom(0))) - HCl1 : wp->atom(2)->distance(*(wp->atom(0))) - HCl2;
                        if(PTC<PTC_cutoff_up and PTC>PTC_cutoff_down){
                            int Hlabel = (HCl1 < HCl2) ? 1 : 2;
                            for(auto const& wp2: wa){
                                if(wp!=wp2){
                                    cout << "test4" << endl ;
                                    double OCl=wp2->atom(0)->distance(Cl);
                                    if ( OCl < OCl_cutoff){
                                        double HCl1 = wp2->atom(1)->distance(Cl);
                                        double HCl2 = wp2->atom(2)->distance(Cl);
                                        //double PTC = (HCl1 < HCl2) ? wp2->atom(1)->distance(*(wp2->atom(0))) - HCl1 : wp2->atom(2)->distance(*(wp2->atom(0))) - HCl2;
                                        int Hlabel2 = (HCl1 < HCl2) ? 1 : 2;
                                        Dp->read(Cl.angle(*(wp->atom(Hlabel)),*(wp2->atom(Hlabel2))));
                                    }
                                }
                            }
                        }
                    }
                }
                cout << "read bead" << fc << " snapshot " << i << '\r' << flush;
            }else{
                cout << "test" << endl ;
                cel.skip(ifs);
                cout << "skip snapshot" << '\r';
            }
        }
    }
    
    ofstream ofs("PTC" + to_string(PTC_cutoff_up) + " to " + to_string(PTC_cutoff_down) + "H-Cl-H.txt");
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

/*
 
 //Cell_ipi_c cel;
 //ifstream ifs1("/Users/jianhangxu/Study/1Cl/Cl_63H2O_pimd/data.pos_0.xyz");
 
 ifstream ifs1("/Users/jianhangxu/Study/1Cl/Cl_63H2O_PBE_vdW_cpmd/cl.cel");
 ifstream ifs2("/Users/jianhangxu/Study/1Cl/Cl_63H2O_PBE_vdW_cpmd/cl.pos");
 Cell_qecp_c cel({1,63,126},{"Cl","O","H"}),cel2({1,63,126},{"Cl","O","H"});
 cel.read(ifs1,ifs2);
 cel2 = cel;
 cel.read(ifs1,ifs2);
 //cel.write_cell(cout) << endl;
 //cel.write_positions(cout) << endl;
 cout << cel2.atom(1).distance(cel2.atom(2)) << endl;
 cout << cel.atom(1).distance(cel.atom(2)) << endl;
 cout << cel2.atom(1).distance(cel.atom(2)) << endl;
 int i,j,k;
 while(cin >> i){
 cout << "num: " << i << " type: " << cel.atom(i).type() << " position: ";
 cel.write_position(cout,i) << endl;
 }
 
 while(cin >> i >> j >> k){
 cout << "position " << i << " : ";
 cel.write_position(cout,i) << " " << cel.atom(i).position() << endl;
 cout << "position " << j << " : ";
 cel.write_position(cout,j) << endl;
 cout << "position " << k << " : ";
 cel.write_position(cout,k) << endl;
 cout << "distance between " << i << " " << j << " : "<< cel.distance(i,j) << " " << cel.atom(i).distance(cel.atom(j)) << endl;
 cout << "angle between " << j << " " << i << " " << k <<  " : " << cel.angle(i,j,k) << " " << cel.atom(i).angle(cel.atom(j),cel.atom(k)) << endl;
 }
//cout << cel.angle(0,10,59) << endl;

 */

/*
 Cell cel,cel2;
 cel.cell_parameters = Eigen::Vector3d(10,10,10).asDiagonal();
 cel.positions.resize(3,3);
 cel.positions << 0,0,110, 53,0,-20 ,13,-16,0;
 cel.is_frac=false;
 
 cel2 = cel;
 
 cout << "cart: " << endl;
 cout << cel.positions << endl;
 
 cel.to_frac();
 cout << "frac: " << endl;
 cout << cel.positions << endl;
 
 cel.to_cart();
 cout << "cart: " << endl;
 cout << cel.positions << endl;
 
 cel.to_frac();
 cout << "frac: " << endl;
 cout << cel.positions << endl;
 
 cout << "cel2: " << endl;
 cout << cel2.positions << endl;
 
 Atom at1(cel.atom(0)), at2(cel2.atom(0));
 
 cout << "angle 0,1,2:";
 cout << cel.atom_position(0) << endl;
 cout << cel.cal_angle_f(cel.positions.row(0),cel.positions.row(1),cel.positions.row(2),cel.cell_parameters) << endl;
 cout << cel.angle(0,1,2) << endl;
 cout << at1.angle(1,2) << endl;
 
 cout << "distance:";
 cout << cel.cal_distance_f(cel.positions.row(0),cel.positions.row(2),cel.cell_parameters) << endl;
 cout << cel.distance(0,2) << endl;
 cout << at1.distance(2) << endl;
 
 cout << "angle 0,1,2:";
 cout << cel2.positions.row(0) << endl;
 cout << cel2.cal_angle_c(cel2.positions.row(0),cel2.positions.row(1),cel2.positions.row(2),cel2.cell_parameters) << endl;
 cout << cel2.angle(0,1,2) << endl;
 cout << at2.angle(1,2) << endl;
 
 cout << "distance:";
 cout << cel2.cal_distance_c(cel2.positions.row(0),cel2.positions.row(2),cel2.cell_parameters) << endl;
 cout << cel2.distance(0,2) << endl;
 cout << at2.distance(2) << endl;
 
 
 cout << "atom_position 0" << endl;
 cout << cel.atom_position(0) << endl;
 cout << at1.position() << endl;
 cel.positions.row(0) << 1,2,3;
 cout << cel.atom_position(0) << endl;
 cout << at1.position() << endl;
 cout << cel.positions << endl;
 cel.positions *= 10;
 cout << cel.positions << endl;
 cout << at1.position() << endl;
 
 cout << "volume: " << cel.volume() << endl;
 cel.set_box(11,11,10);
 cout << "volume: " << cel.volume() << endl;
 
 return 1;
 */
