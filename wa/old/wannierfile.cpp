#include "wannierfile.h"
#include "input.h"

Wannierfile::Wannierfile()
{
    allocate_WC = false;
    allocate_dipole = false;
    sort_WC = false;
    tot_dipole.set(0,0,0);
}

Wannierfile::Wannierfile(const Cellfile &Cel_in): Wannierfile()
{
    cel = &Cel_in;
}

Wannierfile::~Wannierfile()
{
    if(allocate_WC)
        delete [] WCs;
}


void Wannierfile::Routine()
{
    cout << "Using Wannier Center with type: " << INPUT.type << endl;
    switch (INPUT.type) {
        case 1:
            cout << "Calculate average water dipole moment" << endl;
            break;
        default:
            cout << "Wrong type of Wannier" << endl;
            break;
    }
    
    switch (INPUT.type) {
        case 1:
        {
            
        }
            break;
    }
    
}

int Wannierfile::Read_Vdipolefile(ifstream &ifs, Vector3<double> *vdipole_file, int nwater)
{
    cout << "Reading vdipole files!" << endl;
    file_assert(ifs,"Can not find Vdipolefile to Read!",1);
    
    int i=0;
    Vector3<double> vdipole;
    while(ifs.good())
    {
        vdipole_file[i].set(0,0,0);
        for(int j=0; j<nwater; j++)
        {
            ifs >> vdipole.x >> vdipole.y >> vdipole.z;
            vdipole_file[i] = vdipole_file[i] + vdipole;
            ifs.ignore(200,'\n');
        }
        i++;
    }
    return i-1;
}

int Wannierfile::Read_Vdipolefile(ifstream &ifs, Vector3<double> **vdipole_file, int nwater)
{
    cout << "Reading vdipole files!" << endl;
    file_assert(ifs,"Can not find Vdipolefile to Read!",1);
    
    int i=0;
    Vector3<double> vdipole;
    while(ifs.good())
    {
        cout << "snapshot: " << i+1 << endl;
        for(int j=0; j<nwater; j++)
        {
            ifs >> vdipole.x >> vdipole.y >> vdipole.z;
            vdipole_file[j][i] =  vdipole;
            ifs.ignore(200,'\n');
        }
        i++;
    }
    return i-1;
}

int Wannierfile::Read_Vdipolefile(ifstream &ifs, Vector3<double> *vdipole_file, int nwater, int n)
{
    cout << "Reading vdipole files of Water No." << n << endl;
    file_assert(ifs,"Can not find Vdipolefile to Read!",1);
    
    int i=0;
    Vector3<double> vdipole;
    while(ifs.good())
    {
        vdipole_file[i].set(0,0,0);
        for(int j=0; j<nwater; j++)
        {
            ifs >> vdipole.x >> vdipole.y >> vdipole.z;
            ifs.ignore(200,'\n');
            if(j == n)
                vdipole_file[i] = vdipole_file[i] + vdipole;
        }
        i++;
    }
    return i-1;
}

int Wannierfile::Read_Dipolefile(ifstream &ifs,Vector3<double> **dipole_file, int nwater)
{
    cout << "Reading vdipole files!" << endl;
    file_assert(ifs,"Can not find Vdipolefile to Read!",1);
    
    int i=0;
    Vector3<double> dipole;
    string buffer1,buffer2;
    while(ifs.good())
    {
        ifs >> buffer1 >> buffer2;
        cout << "snapshot: " << buffer1 << endl;
        for(int j=0; j<nwater; j++)
        {
            ifs >> dipole.x >> dipole.y >> dipole.z;
            dipole_file[j][i] =  dipole;
            ifs.ignore(200,'\n');
        }
        i++;
    }
    return i-1;
    }

void Wannierfile::Read_Wannier(Cellfile &Cel)
{
    assert(Cel.allocate_atoms);
    assert(Cel.allocate_wan_centers);
    
    if(!allocate_waters)
        Read_water(Cel);
    
    WCs = new Wannier [nwater];
    allocate_WC = true;
    
    time = Cel.time;
    snapshot = Cel.snapshot;
    
    for(int i=0; i< nwater; i++)
    {
        //can add choice here
        //if( INPUT.system=="hydronium" and water[ia].nH != 3) continue;
        WCs[i].numWC=0;
        WCs[i].dipole.set(0,0,0);
        //find the wannier center near the water
        for(int j=0; j< Cel.nband; j++)
        {
            double dO = Cel.atoms[O].pos[i].distance_BC(Cel.wan_centers[j],Cel.celldm)*INPUT.unitconv;
            //read wannier center
            if( dO < INPUT.OW_distance ) // 0.8 Angstroms is a safe cutoff
            {
                WCs[i].index_WC.push_back(j);
                WCs[i].numWC++;
            }
        }
        // O has 6 positive charge H has one WC has two
        if( WCs[i].numWC != 4 )
        {
            cout << "Something went wrong while reading wannier center in " << i << "th water molecule!" <<endl;
            cout << "It has " << WCs[i].numWC << " wannier centers (should be 4)" << endl;
            cout << "Try to change a wannier center lenth cutoff!"<< endl;
            exit(1);
        }
    }
    
}

void Wannierfile::Sort_Wannier(Cellfile &Cel)
{
    assert(allocate_WC);
    
    for(int i=0; i<nwater; i++)
    {
        for(int j=0; j<waters[i].numH; j++)
        {
            double dmin = Cel.atoms[H].pos[waters[i].indexH[j]].distance_BC(Cel.wan_centers[WCs[i].index_WC[j]],Cel.celldm); // temp distance WC[j] and H[j]
            for(int k=j+1; k<WCs[i].numWC; k++)
            {
                double dtmp = Cel.atoms[H].pos[waters[i].indexH[j]].distance_BC(Cel.wan_centers[WCs[i].index_WC[k]],Cel.celldm); // temp distance WC[k] and H[j]
                if(dtmp < dmin)
                {
                    dmin = dtmp;
                    int tmp = WCs[i].index_WC[k];
                    WCs[i].index_WC[k] = WCs[i].index_WC[j];
                    WCs[i].index_WC[j] = tmp;
                }
            }
        }
    }
    sort_WC = true;
}

void Wannierfile::Read_Dipole(Cellfile &Cel, int type)
{
    assert(allocate_WC);
    tot_dipole.set(0,0,0);
    //calculate the dipole moment of the water
    switch (type) {
        case 0:
            //calculation total dipole
            for(int i=0; i< nwater; i++)
            {
                if(waters[i].numH == 2)
                {
                    WCs[i].dipole.set(0,0,0);
                    for(int k=0; k< waters[i].numH; k++)
                    {
                        WCs[i].dipole = WCs[i].dipole + Cel.atoms[H].pos[waters[i].indexH[k]].shortest(Cel.atoms[O].pos[waters[i].indexO],Cel.celldm);
                    }
                    //waters[i].dipole.print2screen();
                    for(int j=0; j<WCs[i].numWC; j++)
                    {
                        WCs[i].dipole = WCs[i].dipole + (Cel.atoms[O].pos[waters[i].indexO].shortest(Cel.wan_centers[WCs[i].index_WC[j]],Cel.celldm))*2.0;
                    }
                    //set unit to eA
                    WCs[i].dipole = WCs[i].dipole*INPUT.unitconv;
                    // set unit 1Debye=0.20819434eA
                    WCs[i].dipole = WCs[i].dipole/0.20819434;
                }
                tot_dipole = tot_dipole + WCs[i].dipole;
            }
            allocate_dipole = true;
            break;
        case 1:
            //calculate dipole contribution of ions
            for(int i=0; i< nwater; i++)
            {
                if(waters[i].numH == 2)
                {
                    WCs[i].dipole.set(0,0,0);
                    for(int k=0; k< waters[i].numH; k++)
                    {
                        WCs[i].dipole = WCs[i].dipole + Cel.atoms[H].pos[waters[i].indexH[k]].shortest(Cel.atoms[O].pos[waters[i].indexO],Cel.celldm);
                    }
                    //set unit to eA
                    WCs[i].dipole = WCs[i].dipole*INPUT.unitconv;
                    // set unit 1Debye=0.20819434eA
                    WCs[i].dipole = WCs[i].dipole/0.20819434;
                }
                tot_dipole = tot_dipole + WCs[i].dipole;
            }
            allocate_dipole = true;
            break;
        case 2:
            //calculate dipole contribution of electrons
            for(int i=0; i< nwater; i++)
            {
                if(waters[i].numH == 2)
                {
                    WCs[i].dipole.set(0,0,0);
                    for(int j=0; j<WCs[i].numWC; j++)
                    {
                        WCs[i].dipole = WCs[i].dipole + (Cel.atoms[O].pos[waters[i].indexO].shortest(Cel.wan_centers[WCs[i].index_WC[j]],Cel.celldm))*2.0;
                    }
                    //set unit to eA
                    WCs[i].dipole = WCs[i].dipole*INPUT.unitconv;
                    // set unit 1Debye=0.20819434eA
                    WCs[i].dipole = WCs[i].dipole/0.20819434;
                }
                tot_dipole = tot_dipole + WCs[i].dipole;
            }
            allocate_dipole = true;
            break;
        case 3:
            //calculate dipole contribution of lone-pair
            if(!sort_WC)
            {
                Sort_Wannier(Cel);
                cout << "Sort wannier" << endl;
            }
            for(int i=0; i< nwater; i++)
            {
                if(waters[i].numH == 2)
                {
                    WCs[i].dipole.set(0,0,0);
                    for(int j=2; j<WCs[i].numWC; j++)
                    {
                        WCs[i].dipole = WCs[i].dipole + (Cel.atoms[O].pos[waters[i].indexO].shortest(Cel.wan_centers[WCs[i].index_WC[j]],Cel.celldm))*2.0;
                    }
                    //set unit to eA
                    WCs[i].dipole = WCs[i].dipole*INPUT.unitconv;
                    // set unit 1Debye=0.20819434eA
                    WCs[i].dipole = WCs[i].dipole/0.20819434;
                }
                tot_dipole = tot_dipole + WCs[i].dipole;
            }
            allocate_dipole = true;
            break;
        case 4:
            //calculate dipole contribution of bond-pair
            if(!sort_WC)
            {
                Sort_Wannier(Cel);
                cout << "Sort wannier" << endl;
            }
            for(int i=0; i< nwater; i++)
            {
                if(waters[i].numH == 2)
                {
                    WCs[i].dipole.set(0,0,0);
                    for(int j=0; j<waters[i].numH; j++)
                    {
                        WCs[i].dipole = WCs[i].dipole + (Cel.atoms[O].pos[waters[i].indexO].shortest(Cel.wan_centers[WCs[i].index_WC[j]],Cel.celldm))*2.0;
                    }
                    //set unit to eA
                    WCs[i].dipole = WCs[i].dipole*INPUT.unitconv;
                    // set unit 1Debye=0.20819434eA
                    WCs[i].dipole = WCs[i].dipole/0.20819434;
                }
                tot_dipole = tot_dipole + WCs[i].dipole;
            }
            allocate_dipole = true;
            break;
        default:
            break;
    }
    
}

void Wannierfile::Read_Vdipole(Wannierfile &WF_s,Wannierfile &WF,ofstream &ofs)
{
    assert(allocate_WC);
    assert(allocate_dipole);
    file_assert(ofs,"Can not write Vdipole to file!",1);
    
    if(WF.allocate_dipole)
    {
        Vector3<double> vp;
        double dt = (time - WF.time);
        // set unit to a.u.time
        dt *= (100000/2.418884326505);
        
        //print vdipole to ofs
        for(int i=0; i<nwater; i++)
        {
            vp = (WCs[i].dipole - WF.WCs[i].dipole)/dt;
            vp.print(ofs);
            //ofs << " " << vp.norm() << endl;
        }
    }
    else
    {
        if(WF_s.allocate_dipole)
        {
            WF.allocate_dipole = true;
            WF.WCs = new Wannier [nwater];
        }
        else
        {
            WF_s.allocate_dipole = true;
            WF_s.WCs = new Wannier [nwater];
        }
    }
    if(WF.allocate_dipole)
    {
        WF.time = WF_s.time;
        for(int i=0; i<nwater; i++)
        {
            WF.WCs[i].dipole = WF_s.WCs[i].dipole;
        }
    }
    //put this to WF_s to store
    WF_s.time = time;
    for(int i=0; i<nwater; i++)
    {
        WF_s.WCs[i].dipole = WCs[i].dipole;
    }
}

void Wannierfile::Read_Vdipole(Wannierfile &WF_s,Wannierfile &WF,ofstream &ofs,Vector3<double> &tot_vdipole)
{
    assert(allocate_WC);
    assert(allocate_dipole);
    file_assert(ofs,"Can not write Vdipole to file!",1);
    
    if(WF.allocate_dipole)
    {
        Vector3<double> vp;
        double dt = (time - WF.time);
        // set unit to a.u.time
        dt *= (100000/2.418884326505);
        
        //print vdipole to ofs
        for(int i=0; i<nwater; i++)
        {
            vp = (WCs[i].dipole - WF.WCs[i].dipole)/dt;
            vp.print(ofs);
            //ofs << " " << vp.norm() << endl;
        }
        tot_vdipole = (tot_dipole - WF.tot_dipole)/dt;
    }
    else
    {
        if(WF_s.allocate_dipole)
        {
            WF.allocate_dipole = true;
            WF.WCs = new Wannier [nwater];
        }
        else
        {
            WF_s.allocate_dipole = true;
            WF_s.WCs = new Wannier [nwater];
        }
    }
    if(WF.allocate_dipole)
    {
        WF.time = WF_s.time;
        for(int i=0; i<nwater; i++)
        {
            WF.WCs[i].dipole = WF_s.WCs[i].dipole;
        }
        WF.tot_dipole = WF_s.tot_dipole;
    }
    //put this to WF_s to store
    WF_s.time = time;
    for(int i=0; i<nwater; i++)
    {
        WF_s.WCs[i].dipole = WCs[i].dipole;
    }
    WF_s.tot_dipole = tot_dipole;
}

void Wannierfile::Print_Dipole(ofstream &ofs)
{
    assert(allocate_dipole);
    file_assert(ofs,"Can not write Dipole to file!",1);
    
    ofs << snapshot << " " << time << endl;
    for(int i=0; i< nwater; i++)
    {
        WCs[i].dipole.print(ofs,0);
        ofs << "  " << WCs[i].dipole.norm() << endl;
    }
}

//add at 2018.4.6 by jianhang
Vector3<double> Wannierfile::posWC(int nO, int nWC) const
{
    if(nWC>=WCs[nO].numWC)
        throw(std::range_error("WANNIERFILE_POSWC: nH >= numH"));
    return cel->wan_centers[Show_indexWC(nO,nWC)];
}



