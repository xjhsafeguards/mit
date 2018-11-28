#include "interface.h"
#include "input.h"

Interface::Interface()
{
    Gs = 0;
    num_z=0;
    allocate_posz = false;
}

Interface::~Interface()
{
    if(allocate_posz)
    {
        delete [] posz;
        delete [] type_mark;
    }
}

void Interface::Read_Gibbs_surface(const Density &dens)
{
    cout << "Finding Gibbs surface!" <<endl;
    
    assert(dens.allocate_dens);
    
    //cout << INPUT.pos_start[0] << " " <<INPUT.pos_start[1] << endl;
    //find the start position
    int i_start[2];
    for(int k=0; k<2 ;k++)
    {
        bool pos_found = false;
        for(int i=0; i<dens.delta; i++)
        {
            if(dens.position[i] >= INPUT.pos_start[k])
            {
                i_start[k] = i;
                pos_found = true;
                break;
            }
        }
        if(!pos_found)
        {
            cerr << "Cannot find start point for Gibbs surface!" << endl;
            exit(1);
        }
    }
    
    int bound[2]; // mark the boundary where is stable
    
    for(int i=0; i<dens.delta; i++)
    {
        if( abs(dens.dens[(i_start[0]+i)%dens.delta]-dens.dens[i_start[0]]) > INPUT.eps )
        {
            bound[0] = (i_start[0]+i)%dens.delta;
            break;
        }
    }
    for(int i=0; i<dens.delta; i++)
    {
        if(i_start[1]-i < 0) i_start[1] += dens.delta;
        if( abs(dens.dens[(i_start[1]-i)%dens.delta]-dens.dens[i_start[1]%dens.delta]) > INPUT.eps )
        {
            bound[1] = (i_start[1]-i)%dens.delta;
            break;
        }
        
    }
    
    double d0[2] = {dens.dens[bound[0]],dens.dens[bound[1]]};
    double n[2] = {0};
    
    while(bound[0] != bound[1])
    {
        if(n[0]> n[1])
        {
            n[1] += abs(dens.dens[bound[1]] - d0[1]);
            bound[1] = bound[1]-1;
            if(bound[1] < 0) bound[1] += dens.delta;
        }
        else
        {
            n[0] += abs(dens.dens[bound[0]] - d0[0]);
            bound[0] = (bound[0]+1)%dens.delta;
        }
    }
    
    Gs = dens.position[bound[0]];
    cout << "Gibbs surface : " << Gs << endl;
}


int Interface::Read_pos_ILI(ifstream &ifs, bool skip)
{
    if(ifs.eof())
    {
        cout << "Read the last pos_ILI file" << endl;
        return 0;
    }
    assert(ifs.good());
    
    allocate_posz = 1;
    num_z = INPUT.atom_num_tot;
    posz = new double [num_z];
    type_mark = new int [INPUT.ntype];
    
    //set the number of type marker
    for(int i=0; i<INPUT.ntype; i++)
    {
        type_mark[i] = 0;
        for(int j=0;j<i;j++)
            type_mark[i] += INPUT.atom_num[j];
    }
    
    string tmpid;
    int tmpnum;
    
    ifs >> snapshot >> time;
    ifs.ignore(200,'\n');
    int count = 0;

    if(skip)
    {
        cout << "Skip Snapshot: " << snapshot << endl;
        Skip_lines(ifs, INPUT.atom_num_tot);
        return 1;
    }
    
    cout << "Read Snapshot: " << snapshot << endl;
    for(int i=0; i<INPUT.ntype; i++)
    {
        for(int j=0; j<INPUT.atom_num[i]; j++)
        {
            ifs >> tmpid >> tmpnum >> posz[count++];
        }
        if(tmpid != INPUT.atom_name[i] || tmpnum != INPUT.atom_num[i])
        {
            cout << "ILI_file: " << INPUT.in_file << " does not match input_file: " << INPUT.geo_file << endl;
            exit(1);
        }
    }

    return 1;
}
