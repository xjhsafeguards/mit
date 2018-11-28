#ifndef INPUT_H
#define INPUT_H

#include "vec3.h"
#include "gfun.h"

class Input
{
    public:
    
    Input();
    ~Input();
    
    string calculation;
    int type; //in organize water 1 to highlight water
             //in density 1 to simple, 2 to Gs, 3 to ili_density
    string ensemble; // npt nvt
    
    string in_file;
    string cel_file;
    string geo_file;
    string wan_file;
    string out_file;
    
    int ss_start;
    int ss_stop; //snapshot information
    int ss_step;
    int ss_n;
    double dt;
    
    int ntype;
    int nband; //total bands of the system
    int *atom_num;
    int atom_num_tot;// total atom numbers
    string *atom_name;
    double *atom_mass;
    
    double T;//temperature of the system
    Vector3<double> Celldm;
    
    // SI for lenth is A
    
    double unitconv;// (input unit/SI unit)
    double OH_distance;//to find water(SI)
    double OO_distance;//to find H_bond(SI)
    double HOO_angle;//to find H_bond(degree)
    double OW_distance; //to find wannier center near O (SI)
    double HW_distance; //to find wannier center near H (SI)
    
    double cutoff;//to record the cutoff of some calculation
    double alpha;// combine with cutoff to smooth functions
    int thicken;// to increase number of data points
    double upper_limit;//to record the range of output
    int delta; //to divide things up
    double eps; //parameter to control precision(SI)
    double displacement;//to save some needed displacement(input unit)
    double pos_start[2]; // to save some start point for calculation(SI)
    int pos_start_n;// to save how many segments of pos_start
    
    int tcf_max; // from 0,0 tom tcf_max,0
    int index[8] = {0}; //for pick up certain 1
    
    public:
    
    
    void Read_input(string inputfile = "INPUT");
    void Print2screen(void); // for test
    
    /* search functions*/
    
    // add by jianhang 2018.3.30
    int Type_index(string type_name) const; // return the value where atom_name[value] = type_name; -1 represent for did not find
    int Atom_num(string type_name) const; // return atom_num[type_name]
    int Atom_num(int type_index) const; // return atom_num[type_name]

    public:
    ifstream ifs_cel;
    streampos cel_top;
    
};

extern Input INPUT;
// for mpi
extern ofstream ofs_log;
extern int NCOR;
extern int RANK;

inline int Input::Type_index(string type_name) const
{
    for(int i=0; i<ntype; i++)
        if(atom_name[i]==type_name)
            return i;
    throw runtime_error("Can not find type");
    return -1;
}

inline int Input::Atom_num(string type_name) const
{
    return atom_num[Type_index(type_name)];
}

inline int Input::Atom_num(int type_index) const
{
    assert(type_index < ntype);
    return atom_num[type_index];
}

#endif
