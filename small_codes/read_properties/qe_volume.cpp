#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cmath>
using namespace std;

int main(int argc, char ** argv){
    string filename("data.cel");
    string outfile("volume.txt");
    double cellmass(-1),unitconv(1);
    bool fastcal=false;
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            if(strncmp(argv[i],"-n",2) == 0){string tmp=argv[++i];filename=tmp;}
            else if(strncmp(argv[i],"-o",2) == 0){string tmp=argv[++i];outfile=tmp;}
            else if(strncmp(argv[i],"-m",2) == 0){string tmp=argv[++i];cellmass=stod(tmp);cout << "Read in mass: " << cellmass << endl;}
            else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    //open input file
    ifstream ifs(filename);
    assert(ifs.good());
    cout<<"Read in file:"<< filename <<endl;
    //open output file
    ofstream ofs(outfile);
    assert(ofs.good());
    ofs << fixed << setprecision(10) << "#          time       snapshot              volume             density" << endl;
    
    string tmp;
    double time,x,y,z,u=pow(unitconv,3),v=0,total_V=0,d=0,total_D=0;
    int count=0;
    ifs >> tmp >> time;
    if(cellmass<=0){
        while(ifs.good()){
            cout << "capture snapshot:" << ++count << "\r" << flush;
            ofs << setw(15) << tmp << setw(15) << time;
            ifs >> x >> tmp >> tmp >> tmp >> y >> tmp >> tmp >> tmp >> z;
            v = x*y*z*u;
            total_V += v;
            ofs << setw(20) << v << endl;
            ifs.ignore(500,'\n');
            ifs >> tmp >> time;
        }
        ofs << "#average volume: " << setw(20) << total_V/count << endl;
    }else{
        cellmass*=10000/6.022140857/5.2917721067/5.2917721067/5.2917721067;
        while(ifs.good()){
            cout << "capture snapshot:" << ++count << "\r" << flush;
            ofs << setw(15) << tmp << setw(15) << time;
            ifs >> x >> tmp >> tmp >> tmp >> y >> tmp >> tmp >> tmp >> z;
            v = x*y*z*u;
            d = cellmass/v;
            total_D += d;
            total_V += v;
            ofs << setw(20) << v << setw(20) << d << endl;
            ifs.ignore(500,'\n');
            ifs >> tmp >> time;
        }
        ofs << "#average volume and density:  " << setw(20) << total_V/count << setw(20) << total_D/count << endl;
    }

    cout << "Job Done" << endl;
}

