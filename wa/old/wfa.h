//Analyzing waterfile
#ifndef WFA_H
#define WFA_H

#include "waterfile.h"
#include "input.h"

class Wfa;
class Wfat;
ostream& Print_rOH(ostream&,const Wfa&);

class Wfa
{
    friend class Wfat;
    friend ostream& Print_rOH(ostream&,const Wfa&);
    
public:
    
    Wfa(const Waterfile& WF_in, const Cellfile& Cel_in);
    Wfa(const Waterfile& WF_in, const Cellfile& Cel_in, const Input& INPUT_in);
    ~Wfa();
    
    /* basic functions */
    void static Routine();
    void set_unitconv(double d){unitconv = d;}
    double Show_time() const {return Cel.time;}
    
    /* rotaional time correlation function */
    void Read_rOH();
    void Read_rOH(int index_O); // read only 1 water
    void Read_rOH_ofp(int index_O); // read only 1 water
    void Read_rOH(int index_O, int index_H); // read certain OH bond
    ostream& Print_rOH(ostream& os) const;
    
    
private:
    /* basic values */
    const Waterfile& WF;
    const Cellfile& Cel;
    double unitconv; // unitconv from Cel to SI unit
    
    /* rotaional time correlation function */
    int n_OH;
    Vector3<double> rOH; //the average vector from O to H in water covalent bond in SI unit;

};

class Wfat
// in order to analyze t related Wf
{
    
public:
    
    Wfat();
    Wfat(const Input& INPUT_in);
    Wfat(const Wfat& WFAT);
    ~Wfat();
    
    /* operator */
    Wfat& operator=(const Wfat& WFAT);
    
    /* basic functions */
    Wfat& operator<<(const Wfa &wfa);

    void static Routine();
    
    /* rotaional time correlation function */
    Wfat& Read_rOH(const Wfa &wfa); // Read in only directions
    void Set_n_rOH_tcf(int n) {n_rOH_tcf = n;}
    void Calculate_rOH_tcf(); // calculate tcf from 0,0 to n_rOH_tcf*n_jump*dt,0 if possible
    void Calculate_rOH_tcf(int delta); // calculate tcf from 0,0 to delta*n_jump*dt,0 if possible
    void Calculate_rOH_tcf_1();
    void Calculate_t2();
    
    ostream& Print_rOH_tcf(ostream &os = cout) const;
    const vector<double>& Show_rOH_tcf() const {return rOH_tcf;}
    Wfat& Sum_rOH_tcf(const Wfat& WFAT);
    Wfat& Div_rOH_tcf(double d);
    
private:
    
    /* basic values */
    int N;  // number of snapshots included start from 1
    double dt; //time interval of snapshot
    double lt; //time of last snapshot
    //vector<Wfa> Wfas; // wfa for different snapshots
    //vector<double> t; // time for corespond wfa
    
    /* basic funcions */
    
    
    
    /* rotaional time correlation function */
    vector<Vector3<double> > rOH; //hold the unit vector
    bool allocate_rOH_tcf;
    int n_rOH_tcf; // number of tcf calculated not including <0*0>, i.e. last index of rOH_tcf
    int n_jump; // dt_tcf = dt * n_jump
    vector<double> rOH_tcf; //<P2(rOH(t)*rOH(0))>/<P2(rOH(0)^2)> and the first one is <P2(rOH(0)^2)>
    double t2; // rotational correlation time (SI ps)
    
};


inline ostream& Print_rOH(ostream& os,const Wfa& AF)
{
    os << AF.rOH;
    return os;
}

inline ostream& Wfa::Print_rOH(ostream& os) const
{
    os << setw(10) << Cel.time << " " << rOH << endl;
    return os;
}
#endif
