//include WC positions and dipoles
#ifndef WANNIERFILE_H
#define WANNIERFILE_H

#include "wannier.h"
#include "waterfile.h"

class Wannierfile: public Waterfile
{
public:
    
    Wannierfile();
    Wannierfile(const Cellfile &Cel_in);
    ~Wannierfile();
    
    static void Routine();
    
    //basic functions for Ir spectra
    static void Routine_Ir();
    static void Routine_Ir_1(ofstream& ofs, ofstream& ofs_v, int read_type = 0); //Read dipole and vdipole from .geo and .wan, write in ofs ,ofs_v
    static void Routine_Ir_1(ofstream& ofs, ofstream& ofs_v, Vector3<double> *vdipole_file, int read_type = 0); //Read dipole and vdipole from .geo and .wan, write in ofs ,ofs_v and *vdipole_file
    
    //basic functions
    static int Read_Vdipolefile(ifstream &ifs, Vector3<double> *vdipole_file, int nwater);
    static int Read_Vdipolefile(ifstream &ifs, Vector3<double> **vdipole_file, int nwater);
    static int Read_Vdipolefile(ifstream &ifs, Vector3<double> *vdipole_file, int nwater,int n);
    static int Read_Dipolefile(ifstream &ifs,Vector3<double> **dipole_file, int nwater);
        //read total_vdipole from ifs to vdipole_file return numbers of ss_n in inputfile, read the nth water only
    
    void Read_Wannier(Cellfile &Cel); // read numWC and index_WC
    void Sort_Wannier(Cellfile &Cel); // Sort the Wannier center with H1 H2 lone-pair
    void Read_Dipole(Cellfile &Cel, int type=0);  // read Dipoles
    //1 = read Dipoles for ion with positive q
    //2 = read Dipoles for WCs
    void Read_Vdipole(Wannierfile &WF_s,Wannierfile &WF,ofstream &ofs);
    void Read_Vdipole(Wannierfile &WF_s,Wannierfile &WF,ofstream &ofs, Vector3<double> &tot_vdipole);
        // WF is the wannierfile for last2 snapshot
        // WF_s is the wannierfile for last snapshot
    void Print_Dipole(ofstream &ofs); // print the dipole in waters
    
    //information stored
    double time;//(SI)
    int snapshot;
    Vector3<double> tot_dipole;//(Debye)
    
    Wannier *WCs;
    bool allocate_WC; // allocate wannier centers in Waters
    bool allocate_dipole; // allocate dipole in Waters
    bool sort_WC;
    
    //add 2018.4.6 by jianhang
protected:
    
    void test_WC() const {   my_throw(allocate_WC,"WANNIERFILE_POSWC: !allocate_WC");}
    void test_dipole() const {   my_throw(allocate_dipole,"WANNIERFILE_SHOWDIPOLE: !allocate_dipole");}
    
public:
    Vector3<double> posWC(int nO, int nWC) const; //return the position of indexth water's indexth wanniercenter in Cellfile unit
    Vector3<double> Show_tot_diple() const {test_dipole(); return tot_dipole;}
    Vector3<double> Show_dipole(int nO) const {
        test_dipole();
        return WCs[nO].dipole;
    }
    
    int Show_indexWC(int nO, int nWC) const {
        test_WC();
        return WCs[nO].index_WC[nWC];
    }
};

void TCF(Vector3<double> *A,double dt,int ss_n,int delta,int tjump=0, double t_cutoff=0, double alpha_cutoff=0); // output to file
int TCF2list(double *result, Vector3<double> *A,double dt,int ss_n,int delta,int tjump=0, double t_cutoff=0, double alpha_cutoff=0,bool output = true); // output **add on** to result, return number of output (last index +1);
int TCF2list2(double *result, Vector3<double> *A,Vector3<double> *B,double dt,int ss_n,int delta,int tjump=0, double t_cutoff=0, double alpha_cutoff=0,bool output = true);
int TCF2list3(double *result, Vector3<double> *A,Vector3<double> *B,Vector3<double> *posi, Vector3<double> *posj, Vector3<double> *celldm, double pos_start, double pos_end, double dt,int ss_n,int delta,int tjump=0, double t_cutoff=0, double alpha_cutoff=0,bool output = true);
int TCF2list4(double **result, Vector3<double> *A, Vector3<double> *B, Vector3<double> *posi, Vector3<double> *posj, Vector3<double> *celldm, double pos_start, double pos_end, int pos_n, double dt,int ss_n,int delta,int tjump,double t_cutoff,double alpha_cutoff,bool output = true);
int TCF2list5(double *result, Vector3<double> *A, Vector3<double> *B, Vector3<double> *C, Vector3<double> *D,Vector3<double> *posi, Vector3<double> *posj, Vector3<double> *celldm, double pos_start, double pos_end, double dt,int ss_n,int delta,int tjump,double t_cutoff,double alpha_cutoff,bool output );
//void TCF(Vector3<double> *A,double dt,int ss_n,int delta,int tjump);
//A     tot_vdipole for each time point
//dt    interval of config in ps
//ss_n  index of last snapshot in A
//delta numbers of output time points
//tjump how many dt between output time points
//result    list to store result

void FT(double *y, double dt, int n2, string out_file="FT.txt");
//void FT(double *y, double dt, int n2,int delta);
//y     value of TCF for time each point
//dt    interval of each y
//n2    index of last y
//delta numbers of output time points should be the same as n2, if y[n2] is not going to 0

//void vUACF(int TCOR,int TRUN,int TJUMP,double dt, int N, string in_file="waters_vdipole.txt");
//TCOR  time steps i.e.max time window
//TRUN  total config number
//TJUMP interval betwwen different initial points
//dt    interval of config in ps
#endif
