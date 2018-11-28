#include "cellfile.h"
#include "input.h"

Cellfile::Cellfile()
{
    snapshot = -1;
    time = -1;
    
    natom = 0;
    ntype = 0;
    allocate_atoms = false;
    allocate_wan_centers = false;
}

Cellfile::~Cellfile()
{
    if(allocate_atoms)
        delete [] atoms;
    if(allocate_wan_centers)
        delete [] wan_centers;
}

int Cellfile::ReadGeometry(ifstream &ifs,bool skip)
// return 0 means to the end of the file
{
    assert(INPUT.ntype>0);
    allocate_atoms = true;
    ntype = INPUT.ntype;
    atoms = new Atoms [ntype];
    
    for(int i=0; i<ntype; i++)
    {
        atoms[i].id = INPUT.atom_name[i];
        atoms[i].na = INPUT.atom_num[i];
        natom += atoms[i].na;
        atoms[i].mass = INPUT.atom_mass[i];
    }
    
    if(ifs.eof())
    {
        cout << "Read the last Geometry file" << endl;
        return 0;
    }
    if(ifs.fail())
    {
        cout << "Can not read Geometry for cellfile" << endl;
        exit(1);
    }
    
    ifs >> snapshot >> time;
    ifs.ignore(200,'\n');
    set_celldm();
    
    if(skip)
    {
        cout << "Skip Snapshot: " << snapshot << endl;
        Skip_lines(ifs, natom);
    }
    else
    {
        cout << "Read Snapshot: " << snapshot << endl;
        for(int i=0; i<ntype; i++) atoms[i].read_geo(ifs);
    }
    
    return 1;
    
}

void Cellfile::ReadInfile(ifstream &ifs)
{
    allocate_atoms = true;
    ntype = INPUT.ntype;
    set_celldm();
    atoms = new Atoms [ntype];
    
    for(int i=0; i<ntype; i++)
    {
        atoms[i].id = INPUT.atom_name[i];
        atoms[i].na = INPUT.atom_num[i];
        natom += atoms[i].na;
        atoms[i].mass = INPUT.atom_mass[i];
    }
    
    streampos original = ifs.tellg();
    
    if(ifs.fail() || ifs.eof())
    {
        cout << "Can not read In_file for cellfile" << endl;
        exit(1);
    }
    
    
    cout << "Read In_file: " << INPUT.in_file << endl;
    for(int i=0; i<ntype; i++)
    {
        ifs.seekg(original);
        atoms[i].read_pos2(ifs);
    }
}

int Cellfile::ReadWC(ifstream &ifs,bool skip)
{
    assert(INPUT.nband>0);
    allocate_wan_centers = true;
    nband = INPUT.nband;
    
    wan_centers = new Vector3<double> [nband];
    
    if(ifs.eof())
    {
        cout << "Read the last Wanniercenter file" << endl;
        return 0;
    }
    if(ifs.fail())
    {
        cout << "Can not read Wanniercenter for cellfile" << endl;
        exit(1);
    }
    
    int tmpss;
    double tmpt;
    
    ifs >> tmpss >> tmpt;
    ifs.ignore(200,'\n');
    
    if(snapshot !=-1 and snapshot != tmpss)
    {
        cout << "geo_file and wan_file are not compatable" << endl;
        exit(1);
    }
    if(skip)
    {
        cout << "Skip Snapshot: " << tmpss << " (wannier)" << endl;
        Skip_lines(ifs, nband);
    }
    else
    {
        cout << "Read Snapshot: " << tmpss << " (wannier)" << endl;
        for(int i=0; i<nband; i++)
            ifs >> wan_centers[i].x >>  wan_centers[i].y >> wan_centers[i].z;
    }
    return 1;
}

void Cellfile::organize_pos()
{
    assert(allocate_atoms);
    
    for(int i=0; i<ntype; i++)
    {
        for(int j=0; j<atoms[i].na; j++)
        {
            atoms[i].pos[j].normal_BC(celldm);
        }
    }
}

void Cellfile::organize_wan()
{
    assert(allocate_wan_centers);
    
    for(int i=0; i<nband; i++)
        wan_centers[i].normal_BC(celldm);
}

void Cellfile::set_celldm()
{
    //change INPUT.Celldm if ensemble is npt
    if(INPUT.ensemble == "npt")
    {
        assert(INPUT.ifs_cel.good());
        int tmps;
        double tmpt;
        
        INPUT.ifs_cel >> tmps >> tmpt;
        if(tmps != snapshot or tmpt != time)
        {
            cerr << "cel_file: " << INPUT.cel_file << "is not compatable with geo_file: " << INPUT.geo_file << endl;
            exit(0);
        }
        INPUT.ifs_cel >> INPUT.Celldm.x >> tmpt >> tmpt;
        INPUT.ifs_cel >> tmpt >> INPUT.Celldm.y >> tmpt;
        INPUT.ifs_cel >> tmpt >> tmpt >> INPUT.Celldm.z;
    }
    celldm = INPUT.Celldm;
}

double Cellfile::Show_cell_volume() const
{
    return celldm.x*celldm.y*celldm.z;
}

