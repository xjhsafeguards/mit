#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <map>
using namespace std;

void read_cel(ifstream& ifs, ofstream& ofs, double unit);
void ofs_label(ofstream& ofs, int sc, double tc);
ofstream& if2of(ifstream& ifs, ofstream& ofs,double unit,int n);
void write_info(ofstream& of_info, int natom, int snap_count, vector<string>& name);


int main(int argc, char** argv){

  string filename("data.xc.xyz");
  string outprefix("data");
  double unitconv(1),time_step(0.48e-3);
  bool fastcal=false;
  if(argc != 1)
  {
      for(int i=1; i<argc; ++i)
      {
        if(strncmp(argv[i],"-n",2) == 0){string tmp=argv[++i];filename=tmp;}
        else if(strncmp(argv[i],"-p",2) == 0){string tmp=argv[++i];outprefix=tmp;}
        else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
        else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];time_step=stod(tmp);}
        else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
        else{cout << "Read in Unknow tag " << argv[i] << endl;}

          //if(strncmp(argv[i],"count",5) == 0 or strncmp(argv[i],"-c",2) == 0)//{INPUT.calculation="count";}
          //else if(strncmp(argv[i],"help",4) == 0 or strncmp(argv[i],"-h",2) == 0)
          //{Help::Routine(INPUT.calculation); return 1;}
          //else if(strncmp(argv[i],"--cal=",6) == 0)
          //{string tmp = argv[i];INPUT.calculation=tmp.substr(6);}
          //else if(strncmp(argv[i],"--type=",7) == 0)
          //{string tmp = argv[i];convstring(tmp.substr(7),INPUT.type);}
          //else if(strncmp(argv[i],"-t",2) == 0)
          //{i++; string tmp = argv[i];convstring(tmp,INPUT.type);}
          //else if(strncmp(argv[i],"--delta=",8) == 0)
          //{string tmp = argv[i];convstring(tmp.substr(8),INPUT.delta);}
          //else if(strncmp(argv[i],"-d",2) == 0)
          //{string tmp = argv[++i];convstring(tmp,INPUT.type);}
      }
  }

  //open input file
  ifstream ifs(filename);
  assert(ifs.good());
  cout<<"Read in file:"<< filename <<endl;
  //open outputfile to read
  ofstream of_pos(outprefix + ".pos");
  ofstream of_cel(outprefix + ".cel");
  ofstream of_info(outprefix + ".info");
  assert(of_pos.good());
  assert(of_cel.good());
  assert(of_info.good());
  of_pos << fixed << setprecision(10);
  of_cel << fixed << setprecision(10);
  //initialize record Data
  int snap_count(0),natom;
  double time_count(0);
  vector<string> name;
  string tmp;
  ifs >> natom;
  ofs_label(of_cel,snap_count,time_count);
  ofs_label(of_pos,snap_count,time_count);
  read_cel(ifs,of_cel,unitconv);
  for(int i=0; i<natom; ++i){
    ifs >> tmp;
    name.push_back(tmp);
    if2of(ifs,of_pos,unitconv,3) << endl;
    ifs.ignore(500,'\n'); // ignore the words after postion in ifs
  }
  ifs.ignore(500,'\n'); // ignore the natom line in ifs
  while( ifs.good()){
    time_count+=time_step;
    ofs_label(of_cel,++snap_count,time_count);
    ofs_label(of_pos,snap_count,time_count);
    read_cel(ifs,of_cel,unitconv);
    for(int i=0; i<natom; ++i){
      ifs >> tmp;
      assert(name[i]==tmp);
      if2of(ifs,of_pos,unitconv,3) << endl;
      ifs.ignore(500,'\n'); // ignore the words after postion in ifs
    }
    cout << "capture snapshot:" << snap_count << "\r" << flush;
    ifs.ignore(500,'\n'); // ignore the natom line in ifs
  }
  write_info(of_info,natom,snap_count,name);


}

//currently only orthronbic
void read_cel(ifstream& ifs, ofstream& ofs, double unit=1){
  string temp;
  for(int i=0;i<3;++i)
    ifs >> temp;
  ofs << setw(15)<< stod(temp)*unit  << setw(15) << 0.0  << setw(15) << 0.0 << endl;
  ifs >> temp;
  ofs << setw(15)<< 0.0  << setw(15) << stod(temp)*unit  << setw(15) << 0.0 << endl;
  ifs >> temp;
  ofs << setw(15)<< 0.0  << setw(15) << 0.0  << setw(15) << stod(temp)*unit << endl;
  ifs.ignore(500,'\n'); // ignore the words after celldm in ifs
}

void ofs_label(ofstream& ofs, int sc, double tc){
  ofs << sc << "  " << tc << endl;
}

ofstream& if2of(ifstream& ifs, ofstream& ofs,double unit=1,int n=1){
    double tmp;
    for(int i=0;i<n;i++){
        ifs >> tmp;
        ofs << tmp*unit << " ";
    }
    return ofs;
}

void write_info(ofstream& of_info, int natom, int snap_count, vector<string>& name){
  of_info << "num of atoms: " << natom << endl;
  of_info << "num of snapshots " << snap_count << endl;
  of_info << "atoms" << endl;
  map<string,int> name_count;
  for(auto it=name.cbegin(); it!=name.cend(); ++it){
    of_info << *it << endl;
    ++name_count[*it];
  }
  for(auto it=name_count.cbegin(); it!=name_count.cend(); ++it){
    of_info << it->first << ": " << it->second  << endl;
  }
}
