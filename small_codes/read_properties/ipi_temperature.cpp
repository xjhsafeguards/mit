#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cmath>
using namespace std;

vector<string> Read_col(istream& is, int col);

int main(int argc, char ** argv){
    string filename("data.out");
    string outfile("temperature.txt");
    int temp_col(4);
    int time_col(2);
    int step_col(1);
    //double unitconv(1);
    //bool fastcal=false;
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            if(strncmp(argv[i],"-n",2) == 0){string tmp=argv[++i];filename=tmp;}
            else if(strncmp(argv[i],"-o",2) == 0){string tmp=argv[++i];outfile=tmp;}
            else if(strncmp(argv[i],"-c",2) == 0){string tmp=argv[++i];temp_col=stoi(tmp);}
            else if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];time_col=stoi(tmp);}
            else if(strncmp(argv[i],"-s",2) == 0){string tmp=argv[++i];step_col=stoi(tmp);}
            //else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            //else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    //open input file
    ifstream ifs(filename);
    assert(ifs.good());
    cout<<"Read in file: "<< filename <<endl;
    cout << "Capture temp in column " << temp_col << endl;
    //open output file
    ofstream ofs(outfile);
    assert(ofs.good());
    ofs << fixed << setprecision(10) << "#       snapshot             time      temperature" << endl;
    
    string tmp;
    auto mark = ifs.tellg();
    
    while(ifs.good()){
        ifs >> tmp;
        if(tmp == "#"){
            ifs.ignore(500,'\n');
            mark = ifs.tellg();
        }
        else
            break;
    }
    
    ifs.seekg(mark);
    vector<string> step= Read_col(ifs,step_col);
    cout << step.size() << endl;
    ifs.clear();
    ifs.seekg(mark);
    vector<string> time= Read_col(ifs,time_col);
    cout << time.size() << endl;
    ifs.clear();
    ifs.seekg(mark);
    vector<string> temp= Read_col(ifs,temp_col);
    cout << temp.size() << endl;
    
    assert(step.size()==time.size());
    assert(step.size()==temp.size());
    
    int count(0);
    double sum(0);
    auto it1 = step.cbegin();
    auto it2 = time.cbegin();
    auto it3 = temp.cbegin();
    for(; it1!=step.cend(); ++it1,++it2,++it3){
        ++count;
        sum+=stod(*it3);
        ofs << setw(15) << static_cast<int>(stod(*it1)) << setw(17) << stod(*it2) << setw(18) << stod(*it3) << '\n';
    }
    ofs << "average temperature: " << sum/count << endl;
    
}

vector<string> Read_col(istream& is, int col){
    vector<string> result;
    string tmp;
    for(int i=0; i<col; ++i)
        is >> tmp;
    while(is.good()){
        result.push_back(tmp);
        is.ignore(5000,'\n');
        for(int i=0; i<col; ++i)
            is >> tmp;
    }
    return std::move(result);
}







