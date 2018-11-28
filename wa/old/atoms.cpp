#include "atoms.h"


Atoms::Atoms()
{
    na=0;
    mass=0;
    allocate_pos=false;
    allocate_index=false;
    pos=NULL;
}

Atoms::~Atoms()
{
    if(allocate_pos) delete [] pos;
    if(allocate_index) delete [] index;
}

void Atoms::read_pos(ifstream &ifs, bool frac)
//ifs start from the the point with typename,x,y,z. Read the same atoms in a line
{
    assert(na>0);
    string idtmp;
    pos = new Vector3<double>[na];
    allocate_pos = true;
    for(int i=0;i<na;i++)
    {
        ifs >> idtmp;
        if (i==0) id = idtmp;
        if (idtmp != id) //quit if read different atoms
        {
            cout << " Reading atom " << idtmp << " is not " << id << endl;
            cerr << "Error in reading atom position";
            exit(1);
        }
        if(frac) //read in fractional coordinate
        {
            ifs >> pos[i].x >> pos[i].y >> pos[i].z;
            
        }
        else // cartesian coordinate
        {
            ifs >> pos[i].x >> pos[i].y >> pos[i].z;
            //cout << " Read " << id << " No." << i+1 << endl;
        }
        if( ifs.fail() )
        {
            cout << " Reading atom " << i+1 << endl;
            cerr << "Error in reading atom position";
            exit(1);
        }
    }
}

void Atoms::read_pos2(ifstream &ifs,bool frac)
//ifs start from the the point with typename,x,y,z. read only certain atomtype of atoms
{
    assert(na>0);
    string idtmp;
    string atomtype = id;
    pos = new Vector3<double>[na];
    index = new int [na];
    allocate_pos = true;
    //allocate_pos = true;
    
    int k = 0;
    for(int i=0;i<na;i++)
    {
        ifs >> idtmp;
        k++;
        while (idtmp != atomtype) //skip if read different atoms
        {
            ifs.ignore(200,'\n');
            k++;
            ifs >> idtmp;
            if(ifs.eof())
            {
                cout << "Can not find enough Atom:" << id << endl;
                exit(1);
            }
        }
        if(frac) //read in fractional coordinate
        {
            ifs >> pos[i].x >> pos[i].y >> pos[i].z;
            index[i] = k;
            ifs.ignore(200,'\n');
            
        }
        else // cartesian coordinate
        {
            ifs >> pos[i].x >> pos[i].y >> pos[i].z;
            index[i] = k;
            ifs.ignore(200,'\n');
            //cout << " Read " << id << " No." << k << "  ";
            //pos[i].print2screen();
        }
        if( ifs.fail() )
        {
            cout << " Reading atom " << id << " No." << i << endl;
            cerr << "Error in reading atom position";
            exit(1);
        }
    }
}

void Atoms::read_geo(ifstream &ifs)
{
    assert(na>0);
    pos = new Vector3<double>[na];
    allocate_pos = true;
    
    for(int i=0;i<na;i++)
    {
        if( ifs.fail() )
        {
            cout << " Reading atom " << i+1 << endl;
            cerr << "Error in reading atom position";
            exit(1);
        }
        
        ifs >> pos[i].x >> pos[i].y >> pos[i].z;
        //cout << pos[i].x << setw(10) << pos[i].y << setw(10) <<pos[i].z << endl;
    }
}
