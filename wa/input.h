#ifndef INPUT_H
#define INPUT_H

#include "gheader.h"
#include "gphy.h"
#include "gfile.h"
#include "ccellp.h"

struct Input;
Inputfile& operator>>(Inputfile&, Input&);

struct Input: public Data{
    
    friend Inputfile& operator>>(Inputfile&, Input&);

    string      calculation     ="test";
    int         type            =0;
    string      ensemble        ="nvt";
    string      system          ="QE";
    
    Inputfile   geo_file;
    Inputfile   cel_file;
    Inputfile   wan_file;

    vector<Inputfile>   in_files;
    vector<File>        out_files;
    string      out_files_type;
    
    //for readcell
    unsigned    ss_start        =1;
    unsigned    ss_stop         =-1;
    unsigned    ss_step         =1;
    
    int         ss_n            =-1;
    double      dt              =-1;
    
    //for cell
    int         ntype           =0;
    int         nband           =0; //total bands of the system
    vector<int>     atom_num;
    vector<string>  atom_name;
    vector<double>  atom_mass;
    vector<double>  atom_charge;
    Cellp           cellp;
    
    //for run
    double      T               =0;//temperature of the system
    //for water and Hbond
    //double      OH_distance     =1.26;//to find water(SI)
    //double      OO_distance     =3.5;//to find H_bond(SI)
    //double      HOO_angle       =30;//to find H_bond(degree)
    //double      OW_distance     =0.8; //to find wannier center near O (SI)
    double      HW_distance     =0.8; //to find wannier center near H (SI)
    
    double      cutoff          =-1;//to record the cutoff of some calculation
    double      alpha           =-1;// combine with cutoff to smooth functions
    int         thicken         =-1;// to increase number of data points
    double      upper_limit     =-1;//to record the range of output
    int         delta           =-1; //to divide things up
    double      eps             =-1; //parameter to control precision(SI)
    double      displacement    =-1;//to save some needed displacement(input unit)
    
    int         tcf_max         =-1; // from 0,0 tom tcf_max,0
    vector<int>     index; //for pick up certain 1
    vector<double>  parameter;
    
    
public:
    
    Input() = default;
    //input
    void read_input(istream& is);
    //output
    int get_atom_num() const;// total atom numbers
    int get_atom_num(size_t i) const;// ith atom numbers
    int get_atom_num(string name) const;
    size_t get_atom_index(string name) const;
    
private:
    //for read
    string unit;
    string buffer;
    bool read;
    //tools
    void test_read(istream& is) const;
    
    template <typename T>
    void read_value(string key, istream& is, T& data){
        if(buffer == key and !read){
            test_read(is);
            getline(is,buffer);
            istringstream iss(buffer);
            iss >> data;
            log() << setw(15) << key << setw(15) << data << endl;
            read = true;
        }
    }
    template <typename T>
    void read_uvalue(string key, istream& is, T& data){
        if(buffer == key and !read){
            test_read(is);
            getline(is,buffer);
            istringstream iss(buffer);
            iss >> data;
            iss >> unit;
            data = data*Unit::unitconv(unit);
            unit = 1;
            log() << setw(15) << key << setw(15) << data << endl;
            read = true;
        }
    }
    template <typename T>
    void read_values(string key, istream& is, vector<T>& data){
        if(buffer == key and !read){
            test_read(is);
            getline(is,buffer);
            istringstream iss(buffer);
            T buffer2;
            log() << setw(15) << key;
            while(iss >> buffer2){
                data.push_back(buffer2);
                log() << setw(15) << buffer2;
            }
            log() << endl;
            read = true;
        }
    }
    
    void read_values(string key, istream& is, vector<Inputfile>& data);
    void read_uvalue(string key, istream& is, Cellp& data);
};

inline void Input::test_read(istream& is) const{
    test(is.good(),"read::INPUT");
}

extern Input INPUT;

#endif
