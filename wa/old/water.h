#ifndef WATER_H
#define WATER_H

#include "gfun.h"
#include "vec3.h"
#include "cellfile.h"

class Water
{
public:
    
    Water();
    ~Water();
    
    int indexO; //index of O, starting from 0
    
    int numH; //total number of H
    int indexH[8]; //index of atom H, starting from 0
    double disH[8]; //distance between H and O
    
public:
    //add jianhang 2018.3.30
    int Show_nH() const {return numH;}
    int Show_indexH(int nH) const {return indexH[nH];}
    
};

//did not consider bondary conditions
class Waterana: public Water
{
public:
    
    Waterana();
    Waterana(const Water &water);
    Waterana(const Water &water,const Cellfile &cel,const double unitconv =1);
    Waterana(const Vector3<double> &inposO,const vector<Vector3<double> > &inposH);
    Waterana(const Waterana &water);
    ~Waterana();
    
    bool allocate_pos;
    Vector3<double> posO;
    vector<Vector3<double> > posH;
    vector<Vector3<double> > OH;
    int set_OH(); //set OH if havent return numbers of OH set
    
    //for 2H water
    Vector3<double> vnorm;
    bool allocate_vnorm;
    int set_vnorm(bool force_set = false); //set vnorm if havent and return 1 or return 0
    
    Waterana translate(Vector3<double> d= Vector3<double>(0.01,0.01,0.01), double mO = 16, double mH = 2);
    Waterana liberate(double d=0.01, double mO = 16, double mH = 2);
    Waterana liberate(Vector3<double> axis, double d=0.01, double mO = 16, double mH = 2);
    Waterana bend(double d=0.01, double mO = 16, double mH = 2);
    Waterana nbend(int steps=1, double d=0.01, double mO = 16, double mH = 2, bool print = false, ostream &os=cout);
    Waterana astretch(double d=0.01, double mO = 16, double mH = 2);
    Waterana sstretch(double d=0.01, double mO = 16, double mH = 2);
    
    //general funtions
    Waterana& operator= (const Waterana &water);
    void set(const Waterana &water);
    void print_xyz(ostream &os = cout,string mark=""); //print current water to os as .xyz format
};

#endif
