#include "rdf.h"
#include "input.h"

Rdf::Rdf()
{
    upper_limit = INPUT.upper_limit;
    delta = INPUT.delta;
    dr = upper_limit/delta;
    allocate_rdf = false;
}

Rdf::~Rdf()
{
    if(allocate_rdf)
        delete [] rdf;
}

void Rdf::Routine()
{
    cout << "Calculating RDF using type: " << INPUT.type << endl;
    switch (INPUT.type) {
        case 0:
            cout << "Reading RDF for OO and OH in 1D" << endl;
            break;
        case 1:
            cout << "Reading RDF for OO and OH in 3D" << endl;
            break;
        case 2:
            cout << "Reading RDF for OMLWF in 1D" << endl;
            break;
        case 3:
            cout << "Reading RDF for OMLWF in 3D" << endl;
            break;
        default:
            cout << "Wrong type of Rdf" << endl;
            exit(1);
            break;
    }
    
    switch (INPUT.type) {
        case 0:
        case 1:
        {
            ifstream ifs(INPUT.geo_file.c_str());
            file_assert(ifs,INPUT.geo_file);
            
            Rdf Ave_OO,Ave_OH;
            Ave_OO.allocate();
            Ave_OH.allocate();
            INPUT.ss_n = 0;
            
            for(int i=1; i<INPUT.ss_stop; i++)
            {
                Cellfile cel;
                
                if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
                {
                    cel.ReadGeometry(ifs,1); // skip the one we dont need
                }
                else
                {
                    if(cel.ReadGeometry(ifs)==0) break; // if to the end of the file skip
                    cel.organize_pos();
                    
                    Rdf rdf_OO,rdf_OH;
                    Waterfile WF;
                    rdf_OO.Read_rdf(cel,WF,"O");
                    Ave_OO += rdf_OO;
                    rdf_OH.Read_rdf(cel,WF,"H");
                    Ave_OH += rdf_OH;
                    INPUT.ss_n++;
                }
            }
            cout << INPUT.ss_n << " snapshots are counted! " << endl;
            Ave_OO /= INPUT.ss_n;
            Ave_OH /= INPUT.ss_n;
            
            Ave_OO.Print_rdf("RDF_OO.txt");
            Ave_OH.Print_rdf("RDF_OH.txt");
        }
            break;
        case 2:
        case 3:
        {
            ifstream ifs(INPUT.geo_file.c_str());
            file_assert(ifs,INPUT.geo_file);
            ifstream ifs_wan(INPUT.wan_file.c_str());
            file_assert(ifs_wan,INPUT.wan_file);
            
            Rdf Ave_OW;
            Ave_OW.allocate();
            INPUT.ss_n = 0;
            
            for(int i=1; i<INPUT.ss_stop; i++)
            {
                Cellfile cel;
                
                if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
                {
                    cel.ReadGeometry(ifs,1); // skip the one we dont need
                    cel.ReadWC(ifs_wan,1);
                }
                else
                {
                    if(cel.ReadGeometry(ifs)==0) break; // if to the end of the file skip
                    if(cel.ReadWC(ifs_wan)==0) break;
                    
                    cel.organize_pos();
                    cel.organize_wan();
                    
                    Rdf rdf_OW;
                    Wannierfile WF;
                    rdf_OW.Read_rdf(cel,WF,INPUT.type-2);
                    Ave_OW += rdf_OW;
                    INPUT.ss_n++;
                }
            }
            cout << INPUT.ss_n << " snapshots are counted! " << endl;
            Ave_OW /= INPUT.ss_n;
            
            Ave_OW.Print_rdf("RDF_OW.txt");;
        }
            break;
        default:
            break;
    }
    
}

void Rdf::Read_rdf(Cellfile &Cel, Waterfile &WF, string atom_id)
{
    assert(Cel.allocate_atoms);
    double volume = Cel.celldm.x * Cel.celldm.y * Cel.celldm.z * pow(INPUT.unitconv,3); // volume for the cell in SI A
    assert(!allocate_rdf);
    if(!WF.allocate_waters)
        WF.Read_water(Cel);
    //allocate then initialize
    allocate();
    //find out which atoms to read
    int id1,id2;
    Vector3<double> *p1,*p2;
    double d; // store the distance of p1,p2 in A
    id1 = WF.O;
    if(atom_id == "O")
    {
        id2 = WF.O;
        p1 = Cel.atoms[id1].pos;
        p2 = Cel.atoms[id2].pos;
        //starts to read rdf
        for(int i=0; i<WF.nwater; i++)
        {
            for(int j=0; j<WF.nwater; j++)
            {
                d = p1[WF.waters[i].indexO].distance_BC(p2[WF.waters[j].indexO],Cel.celldm)*INPUT.unitconv;
                if(d<upper_limit)
                {
                    rdf[(int)(d/dr)] += 1;
                }
                
            }
        }
    }
    if(atom_id == "H")
    {
        id2 = WF.H;
        p1 = Cel.atoms[id1].pos;
        p2 = Cel.atoms[id2].pos;
        //starts to read rdf
        for(int i=0; i<WF.nwater; i++)
        {
            for(int j=0; j<WF.nwater; j++)
                for(int k=0; k<WF.waters[j].numH; k++)
                {
                    d = p1[WF.waters[i].indexO].distance_BC(p2[WF.waters[j].indexH[k]],Cel.celldm)*INPUT.unitconv;
                    if(d<upper_limit)
                        rdf[(int)(d/dr)] += 1;
                }
        }
    }
    *this/= WF.nwater*(Cel.atoms[id2].na/volume); // calculate g(r) by rho*g(r)*dV = dn(r) both volume are in A^3
    for(int i=0; i<delta; i++)
        if(INPUT.type == 0)
            rdf[i]/=dr;
        else if(INPUT.type == 1)
            rdf[i]/= 4*PI*(pow(dr*(i+1),3)-pow(dr*i,3))/3;
}

void Rdf::Read_rdf(Cellfile &Cel,Wannierfile &WF,int outtype)
{
    assert(Cel.allocate_atoms);
    assert(Cel.allocate_wan_centers);
    double volume = Cel.celldm.x * Cel.celldm.y * Cel.celldm.z * pow(INPUT.unitconv,3); // volume for the cell in SI A
    assert(!allocate_rdf);
    if(!WF.allocate_waters)
        WF.Read_water(Cel);
    //allocate then initialize
    allocate();
    //find out which atoms to read
    int id1;
    Vector3<double> *p1,*p2;
    double d; // store the distance of p1,p2 in A
    id1 = WF.O;
    p1 = Cel.atoms[id1].pos;
    p2 = Cel.wan_centers;
    //starts to read rdf
    for(int i=0; i<WF.nwater; i++)
    {
        for(int j=0; j<Cel.nband; j++)
        {
            d = p1[WF.waters[i].indexO].distance_BC(p2[j],Cel.celldm)*INPUT.unitconv;
            if(d<upper_limit)
                rdf[(int)(d/dr)] += 1;
        }
    }
    *this/= WF.nwater; // calculate g(r) by rho*g(r)*dV = dn(r) both volume are in A^3
    if(outtype == 0)
        for(int i=0; i<delta; i++)
            rdf[i]/=dr;
    else if(outtype == 1)
        for(int i=0; i<delta; i++)
            rdf[i]/= (4*PI*(pow(dr*(i+1),3)-pow(dr*i,3))/3)*(Cel.nband/volume);
}

void Rdf::Print_rdf(string out_file)
{
    assert(allocate_rdf);
    ofstream ofs(out_file.c_str());
    file_assert(ofs,out_file);
    
    for(int i=0; i< delta; i++)
        ofs << dr*i << " " << rdf[i] << endl;
}

void Rdf::allocate()
{
    assert(delta!=0);
    rdf = new double [delta];
    allocate_rdf = true;
    for(int i=0; i<delta; i++)
        rdf[i]=0;
}

void Rdf::operator+=(const Rdf &RDF)
{
    assert(delta == RDF.delta);
    assert(upper_limit == RDF.upper_limit);
    for(int i=0; i<delta; i++)
        rdf[i]+=RDF.rdf[i];
}

void Rdf::operator/=(const double x)
{
    assert(allocate_rdf);
    for(int i=0; i<delta; i++)
        rdf[i]/=x;
}
void Rdf::operator*=(const double x)
{
    assert(allocate_rdf);
    for(int i=0; i<delta; i++)
        rdf[i]*=x;
}

