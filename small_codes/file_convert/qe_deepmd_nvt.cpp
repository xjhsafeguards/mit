#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
using namespace std;

void comp_snap(ifstream& ifs, int snapshots);
ofstream& if2of(ifstream& ifs, ofstream& ofs,double unit,int n);

int main(int argc, char** argv){
    
    string filename("water");
    int n_water(63);
    int nO(n_water),nH(n_water*2+1),nforce(n_water*9+3);
    if(argc>1)
        filename = argv[1];
    cout << "set file name to:" << filename << endl;
    ofstream ofs("type.raw");
    for(int i=0;i<nO;i++)
        ofs << "0 ";
    for(int i=0;i<nH;i++)
        ofs << "1 ";
    ofs << endl;
    
    double bs1=23.517,bs2=23.517,bs3=23.517;	
    double force=51.421,energy=27.211385,length=0.52917721092;
  
    //ifstream ifs1(filename+".cel");
    ofstream ofs1("box.raw",std::ofstream::out | std::ofstream::app);
    ifstream ifs2(filename+".evp");
    ofstream ofs2("energy.raw",std::ofstream::out | std::ofstream::app);
    ifstream ifs3(filename+".for");
    ofstream ofs3("force.raw",std::ofstream::out | std::ofstream::app);
    ifstream ifs4(filename+".pos");
    ofstream ofs4("coord.raw",std::ofstream::out | std::ofstream::app);
    //assert(ifs1.good());
    assert(ifs2.good());
    assert(ifs3.good());
    assert(ifs4.good());
    ofs1 << setprecision(10);
    ofs2 << setprecision(10);
    ofs3 << setprecision(10);
    ofs4 << setprecision(10);
    int snapshots,count=0;
    string tmp;
    ifs2 >> snapshots;
    while (ifs2.good() and ifs3.good() and ifs4.good()){
        //comp_snap(ifs2,snapshots);
        comp_snap(ifs3,snapshots);
        comp_snap(ifs4,snapshots);
        //ifs1.ignore(500,'\n');
        //if2of(ifs1,ofs1,length,9) << endl;
        ofs1 << bs1*length << " 0" << " 0" << " 0 " << bs2*length << " 0" << " 0" << " 0 " << bs3*length << endl;
        ifs2 >> tmp;
        ifs2 >> tmp;
        ifs2 >> tmp;
        ifs2 >> tmp;
        ifs2 >> tmp;
        ofs2 << stod(tmp)*energy << endl;
        ifs3.ignore(500,'\n');
        if2of(ifs3,ofs3,force,nforce) << endl;
        ifs4.ignore(500,'\n');
        if2of(ifs4,ofs4,length,nforce) << endl;
        ++count;
        cout << "capture snapshot:" << snapshots << "\r" << flush;
        //ifs1.ignore(500,'\n');
        ifs2.ignore(500,'\n');
        ifs3.ignore(500,'\n');
        ifs4.ignore(500,'\n');
        ifs2 >> snapshots;
    }
    cout << endl;
    cout << count << " snapshots were captured" << endl;
}

void comp_snap(ifstream& ifs, int snapshots){
    int tmp;
    ifs >> tmp;
    assert( tmp == snapshots);
}

ofstream& if2of(ifstream& ifs, ofstream& ofs,double unit=1,int n=1){
    double tmp;
    for(int i=0;i<n;i++){
        ifs >> tmp;
        ofs << tmp*unit << " ";
    }
    return ofs;
}
