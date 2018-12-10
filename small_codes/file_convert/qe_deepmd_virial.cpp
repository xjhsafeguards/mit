#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
using namespace std;

void comp_snap(ifstream& ifs, int snapshots, string filename);
ofstream& if2of(ifstream& ifs, ofstream& ofs,double unit,int n);
ofstream& if2of_cel(ifstream& ifs, ofstream& ofs, double& cel_size, double unit);
int main(int argc, char** argv){
    
    string filename("water");
    int n_water(64);
    int nO(n_water),nH(n_water*2),nforce(n_water*9);
    if(argc>1)
        filename = argv[1];
    cout << "set file name to:" << filename << endl;
    ofstream ofs("type.raw");
    for(int i=0;i<nO;i++)
        ofs << "0 ";
    for(int i=0;i<nH;i++)
        ofs << "1 ";
    ofs << endl;
	
    double force=51.421,energy=27.211385,length=0.52917721092,stress=160.2176621;
    
    ifstream ifs1(filename+".cel");
    ofstream ofs1("box.raw",std::ofstream::out | std::ofstream::app);
    ifstream ifs2(filename+".evp");
    ofstream ofs2("energy.raw",std::ofstream::out | std::ofstream::app);
    ifstream ifs3(filename+".for");
    ofstream ofs3("force.raw",std::ofstream::out | std::ofstream::app);
    ifstream ifs4(filename+".pos");
    ofstream ofs4("coord.raw",std::ofstream::out | std::ofstream::app);
    ifstream ifs5(filename+".str");
    ofstream ofs5("virial.raw",std::ofstream::out | std::ofstream::app);
    assert(ifs1.good());
    assert(ifs2.good());
    assert(ifs3.good());
    assert(ifs4.good());
    assert(ifs5.good());
    ofs1 << setprecision(10);
    ofs2 << setprecision(10);
    ofs3 << setprecision(10);
    ofs4 << setprecision(10);
    ofs5 << setprecision(10);
    int snapshots,count=0;
    string tmp;
    double cel_size;
    ifs1 >> snapshots;
    while (ifs1.good() and ifs2.good() and ifs3.good() and ifs4.good()){
        comp_snap(ifs2,snapshots,"evp");
        comp_snap(ifs3,snapshots,"for");
        comp_snap(ifs4,snapshots,"pos");
        comp_snap(ifs5,snapshots,"str");
        ifs1.ignore(500,'\n');
        if2of_cel(ifs1,ofs1,cel_size,length) << '\n';
        ifs2 >> tmp;
        ifs2 >> tmp;
        ifs2 >> tmp;
        ifs2 >> tmp;
        ifs2 >> tmp;
        ofs2 << stod(tmp)*energy << '\n';
        ifs3.ignore(500,'\n');
        if2of(ifs3,ofs3,force,nforce) << '\n';
        ifs4.ignore(500,'\n');
        if2of(ifs4,ofs4,length,nforce) << '\n';
        ifs5.ignore(500,'\n');
        if2of(ifs5,ofs5,cel_size/stress,9) << '\n';
        ++count;
        cout << "capture snapshot:" << snapshots << "\r"; // << flush;
        ifs1.ignore(500,'\n');
        ifs2.ignore(500,'\n');
        ifs3.ignore(500,'\n');
        ifs4.ignore(500,'\n');
        ifs1 >> snapshots;
    }
    cout << endl;
    cout << count << " snapshots were captured" << endl;
}

void comp_snap(ifstream& ifs, int snapshots, string filename=" "){
    int tmp;
    ifs >> tmp;
    //assert( tmp == snapshots);
    if(tmp != snapshots){
        cerr << filename << " snapshot: " << tmp << " does not match with " << snapshots <<"!" << endl;
        exit(1);
    }
}

ofstream& if2of(ifstream& ifs, ofstream& ofs,double unit=1,int n=1){
    double tmp;
    for(int i=0;i<n;i++){
        ifs >> tmp;
        ofs << tmp*unit << " ";
    }
    return ofs;
}

ofstream& if2of_cel(ifstream& ifs, ofstream& ofs, double& cel_size, double unit=1){
    double tmp;
    cel_size=1;
    ifs >> tmp;
    ofs << tmp*unit << " ";
    cel_size*= tmp*unit;
    ifs >> tmp;
    ofs << tmp*unit << " ";
    ifs >> tmp;
    ofs << tmp*unit << " ";
    ifs >> tmp;
    ofs << tmp*unit << " ";
    ifs >> tmp;
    ofs << tmp*unit << " ";
    cel_size*= tmp*unit;
    ifs >> tmp;
    ofs << tmp*unit << " ";
    ifs >> tmp;
    ofs << tmp*unit << " ";
    ifs >> tmp;
    ofs << tmp*unit << " ";
    ifs >> tmp;
    ofs << tmp*unit << " ";
    cel_size*= tmp*unit;
    return ofs;
}
