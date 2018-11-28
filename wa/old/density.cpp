#include "density.h"
#include "interface.h"
#include "input.h"

Density::Density()
{
    allocate_dens = 0;
    delta = INPUT.delta;
}


Density::~Density()
{
    if(allocate_dens)
    {
        delete [] position;
        delete [] dens;
    }
}

void Density::Routine()
{
    cout << "Compute Density of the System with method: " << INPUT.type << endl;
    switch (INPUT.type) {
        case 0:
            cout << "Calculate both GBS and ILI" << endl;
            break;
        case 1:
            cout << "Calculate only the z axis density distribution" << endl;
            break;
        case 2:
            cout << "Calculate density distribution with GBS" << endl;
            break;
        case 3:
            cout << "Caluclate density distribution with ILI" << endl;
        break;
        default:
        break;
    }
    if(INPUT.type == 0 || INPUT.type == 1 || INPUT.type == 2)
    {
        ifstream ifs(INPUT.geo_file.c_str());
        if(ifs.fail())
        {
            cerr << "Can not open geo_file:" << INPUT.geo_file << endl;
            exit(1);
        }
        streampos original = ifs.tellg(); // save the start point of ifs
        
        ofstream ofs1("density_ss_simple.txt");
        if(ofs1.fail())
        {
            cerr << "Can not open log_file: density_ss_Gs.txt"  << endl;
            exit(1);
        }
        
        Density ave1;
        ave1.allocate();
        ave1.density_simple(ifs,ofs1,1); // only calculate the simplest
        
        if(INPUT.type != 1)                                 // find Gibbs_surface
        {
            Interface interface1;
            interface1.Read_Gibbs_surface(ave1);
            
            // change the format of the output
            for(int i=0;i<ave1.delta;i++)
            {
                if(interface1.Gs < INPUT.pos_start[0])
                {
                    if(ave1.position[i] >= INPUT.pos_start[0])
                        ave1.position[i] -= INPUT.Celldm.z*INPUT.unitconv;
                }
                else
                    if(ave1.position[i] < INPUT.pos_start[0])
                        ave1.position[i] += INPUT.Celldm.z*INPUT.unitconv;
                ave1.position[i] = ave1.position[i] - interface1.Gs;
            }
            ave1.print_density("density_Gs.txt");
        }
    }
    
    if(INPUT.type == 0 || INPUT.type == 3)                  // using ILI
    {
        ifstream ifs1(INPUT.in_file.c_str());
        if(ifs1.fail())
        {
            cerr << "Can not open in_file:" << INPUT.in_file << endl;
            exit(1);
        }
        
        ofstream ofs2("density_ss_ILI.txt");
        if(ofs2.fail())
        {
            cerr << "Can not open log_file: density_ss_ILI.txt"  << endl;
            exit(1);
        }
        
        //ifs.seekg(original);
        Density ave2;
        ave2.allocate(INPUT.displacement,-INPUT.displacement);
        ave2.density_boundary(ifs1,ofs2,INPUT.displacement,1);
    }
    
    
}

void Density::density_boundary(ifstream &ifs,ofstream &ofs,double displacement,bool print)
{
    //Density Ave;
    //Ave.
    //allocate(displacement,-displacement);
    
    int n=0;//count snapshots
    
    for(int i=1; i<INPUT.ss_stop; i++)
    {
        Interface iff;
        Density den;
        den.allocate(displacement,-displacement);
        
        if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
        {
            if(iff.Read_pos_ILI(ifs,1)==0) break; // skip the one we dont need
        }
        else
        {
            if(iff.Read_pos_ILI(ifs)==0) break; // if to the end of the file skip
            iff.organize_posz(displacement, INPUT.Celldm.z);
            den.Read_ILI_simple(iff);
            if(print)
            {
                ofs << "snaptshot: " << iff.snapshot << endl;
                den.print_density(ofs);
            }
            //Ave.
            add(den);
            n++;
        }
        
    }
    
    for(int i=0;i</*Ave.*/delta;i++)
        //Ave.
        dens[i] = /*Ave.*/dens[i]/n;
    if(print)
    {
        //Ave.
        print_density("density_ILI.txt");
    }
    
}


void Density::density_simple(ifstream &ifs, ofstream &logofs,bool print)
{
    //Density Ave; //record Average
    //Ave.
    //allocate();
    
    int n=0;//count snapshots
    
    //cout << Ave.position[0] << endl;
    
    for(int i=1; i<INPUT.ss_stop; i++)
    {
        Cellfile cel;
        Density den;
        den.allocate();
        
        if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
        {
            if(cel.ReadGeometry(ifs,1)==0) break; // skip the one we dont need
        }
        else
        {
            if(cel.ReadGeometry(ifs)==0) break; // if to the end of the file skip
            den.Read_cel_simple(cel);
            if(print)
            {
                logofs << "snaptshot: " << cel.snapshot << endl;
                den.print_density(logofs);           //den.print_density("ABCDEFG");
            }
            
            //Ave.
            add(den);
            n++;
        }
        
    }
    
    for(int i=0;i</*Ave.*/delta;i++)
        //Ave.
        dens[i] = /*Ave.*/dens[i]/n;
    if(print)
        //Ave.
        print_density("density_simple.txt");
    
}


void Density::print_density(string out_file)
{
    ofstream ofs(out_file.c_str());
    if(ofs.fail())
    {
        cerr << "Can not open file: density.txt to write" << endl;
        exit(1);
    }
    for(int i=0;i<delta;i++)
    {
        ofs << position[i] << "\t" << dens[i] << endl;
    }

    ofs.close();
}

void Density::allocate(double pos0,double pos1)
{
    assert(delta>0);
    
    allocate_dens = true;
    position = new double [delta];
    //mass = new double [delta];
    dens = new double [delta];
    delz = (INPUT.Celldm.z - pos0 + pos1)*INPUT.unitconv/delta;
    vol = INPUT.Celldm.x*INPUT.unitconv*INPUT.Celldm.y*INPUT.unitconv*delz;
    
    for(int i=0; i<delta;i++)
    {
        position[i] = pos0*INPUT.unitconv + i*delz;
        dens[i] = 0;
    }
}

void Density::Read_cel_simple(Cellfile &cel)
{
    assert(allocate_dens);
    assert(cel.allocate_atoms);
    
    cel.organize_pos();
    for(int j=0;j<cel.ntype;j++)
    {
        for(int k=0;k<cel.atoms[j].na;k++)
        {
            dens[(int)(cel.atoms[j].pos[k].z*INPUT.unitconv/delz)] += cel.atoms[j].mass/vol/0.6022;
        }
    }
}

void Density::Read_ILI_simple(Interface &iff)
{
    assert(allocate_dens);
    assert(iff.allocate_posz);
    
    for(int i=0; i<INPUT.ntype; i++)
    {
        for(int j=0; j<INPUT.atom_num[i]; j++)
        {
            dens[(int)((iff.get_posz(i,j)*INPUT.unitconv-position[0])/delz)] += INPUT.atom_mass[i]/vol/0.6022;
        }
    }
}
