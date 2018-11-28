#ifndef RDF_H
#define RDF_H

#include "waterfile.h"
#include "wannierfile.h"

class Rdf
{
public:
    Rdf();
    ~Rdf();
    //basic functions
    static void Routine();
    void Read_rdf(Cellfile &Cel,Waterfile &WF,string atom_id); // read rdf for atom O and atom_id
    void Read_rdf(Cellfile &Cel,Wannierfile &WF,int outtype=1); //read rdf for MLWF
    void Print_rdf(string out_file="RDF.txt");
    void allocate();//allocate rdf
    void operator+=(const Rdf &RDF);// add RDF to this
    void operator/=(const double x);
    void operator*=(const double x); // divide each rdf by x

    double *rdf; // cout the average for radial distribution
    bool allocate_rdf;
    //basic parameters
    double upper_limit; // maximum length considered (SI)
    int delta;
private:
    double dr; // step per each delta (SI)
    
};

#endif
