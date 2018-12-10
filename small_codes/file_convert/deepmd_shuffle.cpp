#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <ctime>
#include <string>
#include <fstream>
#include <algorithm>    // std::random_shuffle
#include <cstdlib>      // std::rand, std::srand

using namespace std;

void Read_file(ifstream& ifs,vector<string>& vec){
    assert(ifs.good());
    string tmp;
    getline(ifs,tmp);
    while(ifs.good()){
        vec.push_back(std::move(tmp));
        getline(ifs,tmp);
    }
}

void Cout_progress(int current,int total,int size=20){
    cout << "[";
    int i=0;
    for(;i<current*size/total;++i) cout << "#";
    for(;i<size;++i) cout << " ";
    cout << "] " << setw(5) << current*100/total << "%" << '\r'; 
}

int main(int argc, char** argv){
    string infile_prefix("old_"),outfile_prefix("");
    ifstream ifs1(infile_prefix+"box.raw");
    ifstream ifs2(infile_prefix+"coord.raw");
    ifstream ifs3(infile_prefix+"energy.raw");
    ifstream ifs4(infile_prefix+"force.raw");
    ifstream ifs5(infile_prefix+"virial.raw");
    
    vector<string> box,coord,energy,force,virial;

    Read_file(ifs1,box);
    Read_file(ifs2,coord);
    Read_file(ifs3,energy);
    Read_file(ifs4,force);
    Read_file(ifs5,virial);

    int vsize=box.size();
    assert(vsize == coord.size());
    assert(vsize == energy.size());
    assert(vsize == force.size());
    assert(vsize == virial.size());

    vector<int> rand;
    for(int i=0; i<vsize; ++i) rand.push_back(i);
    random_shuffle(rand.begin(),rand.end());

    ofstream ofs("index.txt");
    ofstream ofs1(outfile_prefix+"box.raw");
    ofstream ofs2(outfile_prefix+"coord.raw");
    ofstream ofs3(outfile_prefix+"energy.raw");
    ofstream ofs4(outfile_prefix+"force.raw");
    ofstream ofs5(outfile_prefix+"virial.raw");
    ofs1.precision(10);
    ofs2.precision(10);
    ofs3.precision(10);
    ofs4.precision(10);
    ofs5.precision(10);
    
    int counting=0;
    for( const auto& i : rand){
        ofs << i << '\n';
        ofs1 << box[i] << '\n';
        ofs2 << coord[i] << '\n';
        ofs3 << energy[i] << '\n';
        ofs4 << force[i] << '\n';
        ofs5 << virial[i] << '\n';
        Cout_progress(++counting,vsize);
    }
    
}
 

