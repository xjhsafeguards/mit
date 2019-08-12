#include <iomanip>
#include <chrono>

#include "Math_all.h"
#include "Cell.h"
#include "Cell_Wannier90.h"
#include "Cell_TEMP.h"
#include "Molecule_water.h"
#include "Utility_time.h"

using namespace std;

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
    double box_size = 12.419790034; //23.47 bohr

#ifdef QUANT
    water_parameter::OH_distance=2.0;
#endif
    Help_content += "-p for set name of position files\n";
    //Help_content += "-fs for the fisrt snapshot to read\n";
    Help_content += "-fe for the last snapshot to read\n";
    //Help_content += "-st for the step of snapshot to read\n";
    Help_content += "-w for set water searching parameter OH_distance\n";
    //Help_content += "-ocl for set the OCl_cutoff\n";
    //Help_content += "-ohl for set the HCl_cutoff\n";
    //Help_content += "-coo for set the OO_cutoff\n";
    Help_content += "-bs for set the set box size\n";
    
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
            else if(strncmp(argv[i],"-bs",3) == 0){string tmp=argv[++i];box_size=stod(tmp);}
            else if(strncmp(argv[i],"-h",2) == 0){cout << Help_content << endl; return 1;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    cout << "Reading folder: " + file_folder << endl;
    cout << "Reading snapshot from: " + to_string(f_start) + " to: " + to_string(f_end) + " with step: " + to_string(f_step) << endl;
    cout << "Box size : " + to_string(box_size) << endl;
    //cout << "OCl_cutoff :" + to_string(OCl_cutoff) << endl;
    //cout << "HCl_cutoff :" + to_string(HCl_cutoff) << endl;
    //cout << "OO_cutoff :" + to_string(OO_cutoff) << endl;
    //cout << "Distribution functions steps: " + to_string(DF_step) + " with gassian sigma(step): " + to_string(DF_sigma) << endl;
    
    ifstream ifs(file_folder);
    
    cell* cel = new cell_temp1();
    cel->set_box(box_size,box_size,box_size);
    molecule_manip* mol = new water_manip();
    
    ofstream ofs("pos_info.txt");
    ofs << setprecision(15);
    ofs << setw(20) << "OH - ClH" << setw(20) << "O-H-Cl angle" << setw(20) << "O-Cl" << endl;
    
    for(int i=0;i!=f_end;++i){
        
        cout << "Reading snapshot: " << i << '\r' << endl;
        cel->clear();
        GTIMER.start("Read Cell");
        if(!ifs.good() or ifs.eof())
            break;
        cel->read_atoms(ifs);
        GTIMER.stop("Read Cell");
        
        
        GTIMER.start("Read Water");
        mol->read(*cel);
        GTIMER.stop("Read Water");

        for(const auto& mol : cel->mols("H2O")){
            double HCl1 = mol->atoms()[1]->distance(*(cel->atoms()[0]));
            double HCl2 = mol->atoms()[2]->distance(*(cel->atoms()[0]));
            if( HCl1 < HCl2 )
                ofs << setw(20) << mol->atoms()[1]->distance(*(mol->atoms()[0])) - HCl1 << setw(20) << mol->atoms()[1]->angle(*(mol->atoms()[0]),*(cel->atoms()[0])) << setw(20) << mol->atoms()[0]->distance(*(cel->atoms()[0])) << setw(10) << mol->atoms().size()-1 << endl;
            else
                ofs << setw(20) << mol->atoms()[2]->distance(*(mol->atoms()[0])) - HCl2 << setw(20) << mol->atoms()[2]->angle(*(mol->atoms()[0]),*(cel->atoms()[0])) << setw(20) << mol->atoms()[0]->distance(*(cel->atoms()[0])) << setw(10) << mol->atoms().size()-1 << endl;
        }
    }
    GTIMER.stop("all");
    GTIMER.summerize(cout);
}
