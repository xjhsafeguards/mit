#include <iomanip>
#include <chrono>

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
    
    string infoldername,outfoldername;
    
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            //if(strncmp(argv[i],"-n",2) == 0){string tmp=argv[++i];filename=tmp;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];snapshot_count=stoi(tmp);}
            if(strncmp(argv[i],"-w",2) == 0){string tmp=argv[++i];water_parameter::OH_distance=stod(tmp);}
            else if(strncmp(argv[i],"-if",3) == 0){string tmp=argv[++i];infoldername=tmp;}
            else if(strncmp(argv[i],"-of",3) == 0){string tmp=argv[++i];outfoldername=tmp;}
            //else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            //else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];time_step=stod(tmp);}
            //else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    ifstream ifs1(infoldername+"cl.cel");
    ifstream ifs2(infoldername+"cl.pos");

    
    for(int i=0;i!=21220;++i){
    //for(int i=0;i!=5011;++i){
        
        std::shared_ptr<cell> cel = make_shared<cell_qecp>(cell_qecp({1,63,126},{"Cl","O","H"}));
        cel->read_box(ifs1);
        cel->read_atoms(ifs2);
        
        if(i>5000 and i%5==0){

            ofstream ofs(outfoldername+to_string(i)+"wannier.tail");
            ofs << setprecision(15);
            ofs << cel->boxp()->data() << endl;
            ofs << "ATOMIC_POSITIONS {angstrom}" << endl;
            for( const auto& at: cel->atoms()){
                ofs << at->get_type() << " " << at->cart() << endl;
            }
            
            
        }
        //     cellv.push_back(cel);
        //cout << "read " << cellv.size() << '\r' << flush;
        cout << "Read snapshot " << i << '\r' << flush;
    }
    
}
    
