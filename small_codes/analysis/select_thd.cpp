#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>

using namespace std;

vector<double> Read_col(istream& is, int col);
double average(const vector<double>& v);
double average( typename vector<double>::const_iterator beg, typename vector<double>::const_iterator end);

int main(int argc, char ** argv){
    string temperature("temperature.txt");
    string hbond("average_Hbs.txt");
    string density("volume.txt");
    string outfile("std.txt");
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            //if(strncmp(argv[i],"-n",2) == 0){string tmp=argv[++i];filename=tmp;}
            //else
            if(strncmp(argv[i],"-t",2) == 0){string tmp=argv[++i];temperature=tmp;}
            else if(strncmp(argv[i],"-h",2) == 0){string tmp=argv[++i];hbond=tmp;}
            else if(strncmp(argv[i],"-d",2) == 0){string tmp=argv[++i];density=tmp;}
            else if(strncmp(argv[i],"-o",2) == 0){string tmp=argv[++i];outfile=tmp;}
            //else if(strncmp(argv[i],"-m",2) == 0){string tmp=argv[++i];cellmass=stod(tmp);cout << "Read in mass: " << cellmass << endl;}
            //else if(strncmp(argv[i],"-angs",5) == 0){assert(!fastcal);unitconv=1.8897161646320723;}
            //else if(strncmp(argv[i],"-f",2) == 0){fastcal=true;}
            else{cout << "Read in Unknow tag " << argv[i] << endl;}
        }
    }
    
    //open input file
    ifstream ifs_temp(temperature);
    assert(ifs_temp.good());
    cout<<"Read in temperature file:"<< temperature <<endl;
    ifstream ifs_hb(hbond);
    assert(ifs_hb.good());
    cout<<"Read in Hbond file:"<< hbond <<endl;
    ifstream ifs_density(density);
    assert(ifs_density.good());
    cout<<"Read in density file:"<< density <<endl;
    
    cout<<"Please make sure that all in files are compatable! " <<endl;
    //open output file
    ofstream ofs(outfile);
    assert(ofs.good());
    ofs << "# calculate ((mean - x)/mean)^2 and sqrt(sum/3)*100" << endl;
    ofs << fixed << setprecision(10) << "#" << setw(9) << "snapshot" << setw(17) << "temperature" << setw(17) << "hb" << setw(17) << "density" << setw(17) << "All" << endl;
    
    vector<double> vtemp(Read_col(ifs_temp,3));
    cout << vtemp.size() << endl;
    vector<double> vhb(Read_col(ifs_hb,5));
    cout << vhb.size() << endl;
    vector<double> vden(Read_col(ifs_density,4));
    cout << vden.size() << endl;
    
    double a1(average(vtemp));
    double a2(average(vhb));
    double a3(average(vden));
    cout << a1 << " " << a2 << " "  << a3 << " " << endl;
    auto it1=vtemp.cbegin();
    auto it2=vhb.cbegin();
    auto it3=vden.cbegin();
    
    int count(0);
    double d1, d2, d3, da;
    while(it1!=vtemp.cend() and it2!=vhb.cend() and it3!=vden.cend()){
        
        d1 = pow((a1-(*it1))/a1,2);
        d2 = pow((a2-(*it2))/a2,2);
        d3 = pow((a3-(*it3))/a3,2);
        da = sqrt((d1 + d2 + d3)/3)*100;
        ofs << setw(10) << ++count << setw(17) << d1 << setw(17) << d2 << setw(17) << d3 << setw(17) << da << '\n';
        ++it1;
        ++it2;
        ++it3;
    }
    
}

vector<double> Read_col(istream& is, int col){
    vector<double> result;
    string tmp;
    stringstream ss;
    string s;
    getline(is,s);
    while(is.good()){
        //cout << s << endl;
        ss.clear();
        ss.str(s);
        for(int i=0; i<col; ++i)
            ss >> tmp;
        //cout << tmp << endl;
        result.push_back(stod(tmp));
        getline(is,s);
    }
    return std::move(result);
}

double average(const vector<double>& v){
    return average(v.cbegin(),v.cend());
}
double average( typename vector<double>::const_iterator beg, typename vector<double>::const_iterator end){
    int count(0);
    double sum(0);
    while(beg!=end){
        sum+=*beg;
        ++count;
        ++beg;
    }
    return sum/count;
}
