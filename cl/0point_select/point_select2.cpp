#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

class structure{
    
    friend ostream& operator<<(ostream&,structure&);
    friend istream& operator>>(istream&,structure&);
    
    int natom = 190;
    double timestep;
    vector<string> lines;
    
public:
    int count_printline(){
        int i=0;
        for(const auto& data : lines){
            stringstream ss(data);
            string tmp;
            vector<string> sd;
            ss >> tmp;
            while(ss.good()){
                sd.push_back(tmp);
                ss >> tmp;
            }
            if(sd.size()==5 and stod(sd[4])<4)
                ++i;
        }
        return i;
    }
};

ostream& operator<<(ostream& os,structure& st){
    for(const auto& data : st.lines){
        /*
        stringstream ss(data);
        string tmp;
        vector<string> sd;
        ss >> tmp;
        while(ss.good()){
            sd.push_back(tmp);
            ss >> tmp;
        }
        if(sd.size()==5 and stod(sd[4])<4)
         */
            os << data << '\n';
    }
    return os;
}
istream& operator>>(istream& is,structure& st){
    string tmp;
    is >> tmp >> st.timestep;
    is.ignore(500,'\n');
    for(int i=0;i<st.natom;++i){
        getline(is,tmp);
        st.lines.push_back(tmp);
    }
    return is;
}

class atom{
    
public:
    double distance;
    double alpha_L,alpha_T1,alpha_T2;
    structure* position;
};

ostream& operator<<(ostream& os,atom& at){
    os << setw(20) << at.distance;
    os << setw(20) << at.alpha_L << setw(20) << at.alpha_T1 << setw(20) << at.alpha_T2;
    return os;
}
istream& operator>>(istream& is,atom& at){
    is >> at.distance;
    is >> at.alpha_L >> at.alpha_T1 >> at.alpha_T2;
    is.ignore(500,'\n');
    return is;
}

int main(){
    
    ifstream ifs1("traj.d/cl.MD.snapshots_1219");
    ifstream ifs2("traj.d/cl.NQEMD");
    
    vector<structure*> data1,data2;
    while(ifs1.good()){
        structure* tmp = new structure();
        ifs1 >> *tmp;
        data1.push_back(tmp);
        //cout << "read 1" << endl;
    }
    while(ifs2.good()){
        structure* tmp = new structure();
        ifs2 >> *tmp;
        data2.push_back(tmp);
        //cout << "read 1" << endl;
    }
    
    ifstream ifs3("alpha.d/cls.eig.dat");
    ifstream ifs4("alpha.d/pi.8beads.eig.dat");
    
    ifs3.ignore(500,'\n');
    ifs4.ignore(500,'\n');
    
    vector<atom*> atoms1,atoms2;
    int i=0;
    while(ifs3.good()){
        atom* tmp = new atom();
        ifs3 >> *tmp;
        tmp->position = data1[i++/64];
        atoms1.push_back(tmp);
    }
    i=0;
    while(ifs4.good()){
        atom* tmp = new atom();
        ifs4 >> *tmp;
        tmp->position = data2[i++/64];
        atoms2.push_back(tmp);
    }
    
    ofstream ofs1("selected_c.xyz");
    ofstream ofs2("selected_q.xyz");
    ofs1 << setprecision(15);
    ofs2 << setprecision(15);
    
    for(const auto& at : atoms1){
        if(abs(at->distance - 3.5) < 0.1 and at->alpha_L < 0.2){
            ofs1 << "190" << '\n';
            //ofs1 << at->position->count_printline() << '\n';
            ofs1 << "# selected classical atom " << *at << '\n';
            ofs1 << *(at->position);
        }
    }
    for(const auto& at : atoms1){
        if(abs(at->distance - 3.5) < 0.1 and at->alpha_L > 0.3){
            ofs1 << "190" << '\n';
            //ofs1 << at->position->count_printline() << '\n';
            ofs1 << "# selected classical atom " << *at << '\n';
            ofs1 << *(at->position);
        }
    }
    for(const auto& at : atoms2){
        if(abs(at->distance - 3.5) < 0.1 and at->alpha_L < 0){
            ofs2 << "190" << '\n';
            //ofs2 << at->position->count_printline() << '\n';
            ofs2 << "# selected NQE atom " << *at << '\n';
            ofs2 << *(at->position);
        }
    }
    for(const auto& at : atoms2){
        if(abs(at->distance - 3.0) < 2 and at->alpha_L > 0.69){
            ofs2 << "190" << '\n';
            //ofs2 << at->position->count_printline() << '\n';
            ofs2 << "# selected NQE atom " << *at << '\n';
            ofs2 << *(at->position);
        }
    }
    
}
