#include "H_bondfile.h"
#include "input.h"

Hbondfile::Hbondfile()
{
    allocate_Hbonds = false;
    avg_aion = 0;
    avg_dion = 0;
    avg_count_ion = 0;
    avg_awater = 0;
    avg_dwater = 0;
    avg_count_water = 0;
}

Hbondfile::~Hbondfile()
{
    if(allocate_Hbonds)
        delete [] Hbonds;
}

void Hbondfile::Routine()
{
    cout << "Analyzing H_bond using type: " << INPUT.type << endl;
    ifstream ifs(INPUT.geo_file.c_str());
    ofstream ofs1("average_Hbs.txt");
    ofstream ofs2("each_Hbs.txt");
    ofstream ofs3("ion_Hbs.txt");
    ofstream ofs4("Hbonds.txt");
    if(ifs.fail())
    {
        cerr << "Can not open in_file:" << INPUT.geo_file << endl;
        exit(1);
    }
    if(ofs1.fail())
    {
        cerr << "Can not open file average_Hbs.txt to write" << endl;
        exit(1);
    }
    ofs1 << "snapshot\ttime\tavg_awater\tavg_dwater\tavg_count_water\tavg_aion\tavg_dion\tavg_count_ion" << endl;
    if(ofs2.fail())
    {
        cerr << "Can not open file each_Hbs.txt to write" << endl;
        exit(1);
    }
    ofs2 << "O\tAccept\tDonate\tA1\tA2\tA3\tA4\tA5\tA6\tA7\tA8\tD1\tD2\tD3\tD4" << endl;
    if(ofs3.fail())
    {
        cerr << "Can not open file ion_Hbs.txt to write" << endl;
        exit(1);
    }
    ofs3 << "ss\ttime\tO\tnumH\tAccept\tDonate\tA1\tA2\tA3\tA4\tA5\tA6\tA7\tA8\tD1\tD2\tD3\tD4" << endl;
    if(ofs4.fail())
    {
        cerr << "Can not open file Hbonds.txt to write" << endl;
        exit(1);
    }
    
    //Hbondfile avg;
    int count=0;
    int countion=0;
    int count_iona[10]={0};
    int count_iond[10]={0};
    int count_water[10]={0};
    //Calculate for each snapshot
    for(int i=1; i<INPUT.ss_stop; i++)
    {
        Cellfile cel;
        if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
            cel.ReadGeometry(ifs,1); // skip the one we dont need
        else
        {
            if(cel.ReadGeometry(ifs)==0) break; // if to the end of the file skip
            cel.organize_pos();
            //generate Hbondfile
            Hbondfile HF;
            HF.Read_Hbond(cel);
            count+=HF.nwater;
            //output average hbs information for one snapshot
            ofs1 << cel.snapshot << '\t' << cel.time << '\t';
            HF.Print_average(ofs1);
            //output hbs information for each water molecule
            ofs2 << cel.snapshot << '\t' << cel.time << endl;
            HF.Print_each(ofs2);
            //collecting overall information
            for(int j=0; j<HF.nwater; j++)
                count_water[HF.Hbonds[j].HB_count]++;
            //output the ion information
            if(INPUT.type == 0)
            {
                for(int j=0; j<HF.nion; j++)
                {
                    ofs3 << cel.snapshot << '\t' << cel.time << '\t' << HF.ions[j] << '\t' << HF.waters[HF.ions[j]].numH << '\t';
                    HF.Hbonds[HF.ions[j]].Print(ofs3);
                    count_water[HF.Hbonds[HF.ions[j]].HB_count]--;
                    count_iona[HF.Hbonds[HF.ions[j]].accept_count]++;
                    count_iond[HF.Hbonds[HF.ions[j]].donate_count]++;
                    countion++;
                    count--;
                }
            }
        }
    }
    double totaion=0;
    double totdion=0;
    double tothwater=0;
    if(INPUT.type == 0)
    {
        ofs4 << " --- THE H-BOND NETWORK OF IONS --- " << endl;
        ofs4 << setw(10) << " HB_number" << setw(20) << "accept_count" << setw(20) << "percentage(%)" << endl;
        for(int i=0;i<10;i++)
        {
            ofs4 << setw(10) << i << setw(20) << count_iona[i] << setw(20) << ((double)count_iona[i])/countion*100 << endl;
            totaion += i*((double)count_iona[i]);
        }
        ofs4 << setw(10) << " HB_number" << setw(20) << "donate_count" << setw(20) << "percentage(%)" << endl;
        for(int i=0;i<10;i++)
        {
            ofs4 << setw(10) << i << setw(20) << count_iond[i] << setw(20) << ((double)count_iond[i])/countion*100 << endl;
            totdion += i*((double)count_iond[i]);
        }
        ofs4 << " Counting for ions = " << countion << endl;
        ofs4 << " Average accept HBs for ions is " << totaion/countion << endl;
        ofs4 << " Average donate HBs for ions is " << totdion/countion << endl;
    }
    ofs4 << " --- THE H-BOND NETWORK OF WATER MOLECLUES --- " << endl;
    ofs4 << setw(10) << " HB_number" << setw(20) << "HB_count" << setw(20) <<"percentage(%)" << endl;
    for(int i=0;i<10;i++)
    {
        ofs4 << setw(10) << i << setw(20) << count_water[i] << setw(20) << ((double)count_water[i])/count*100 << endl;
        tothwater += i*((double)count_water[i]);
    }
    ofs4 << " Average accept HBs for water molecules is " << (tothwater-totaion+totdion)/2/count << endl;
    ofs4 << " Average donate HBs for water molecules is " << (tothwater+totaion-totdion)/2/count << endl;
}

void Hbondfile::Read_Hbond(Cellfile &Cel)
{
    assert(Cel.allocate_atoms);
    
    if(!allocate_waters)
        Read_water(Cel);
    
    Hbonds = new Hbond [nwater];
    allocate_Hbonds = true;
    
    //for i!=j check all the h_bond informations
    for(int i=0; i<nwater; i++)
    {
        for(int j=0; j<i; j++)
        {
            Donate(Cel,i,j);
        }
        for(int j=i+1; j<nwater; j++)
        {
            Donate(Cel,i,j);
        }
    }
    if(nion != 0)
    {
        avg_aion = avg_aion/nion;
        avg_dion = avg_dion/nion;
        avg_count_ion = avg_dion + avg_aion;
    }
    else
    {
        avg_aion = -1;
        avg_dion = -1;
        avg_count_ion = -1;
    }
    avg_awater = avg_awater/(nwater-nion);
    avg_dwater = avg_dwater/(nwater-nion);
    avg_count_water = avg_awater + avg_dwater;
    
}

void Hbondfile::Print_average(ofstream &ofs)
{
    assert(ofs.good());
    ofs << avg_awater << '\t' << avg_dwater << '\t' << avg_count_water << '\t' << avg_aion << '\t' << avg_dion << '\t' << avg_count_ion << '\t' << endl;
}

void Hbondfile::Print_each(ofstream &ofs)
{
    assert(ofs.good());
    assert(allocate_Hbonds);
    
    for(int i=0; i<nwater; i++)
    {
        ofs << i << '\t';
        //print each water
        Hbonds[i].Print(ofs);
    }
}


void Hbondfile::Donate(Cellfile &Cel,int O1, int O2)
{
    assert(allocate_waters);
    assert(allocate_Hbonds);
    assert(O1!=O2);
    
    //shortcut for the postion of O and H used next
    Vector3<double> pO1,pO2,pH;
    pO1 = Cel.atoms[O].pos[O1];
    pO2 = Cel.atoms[O].pos[O2];
    
    if(pO1.distance_BC(pO2,Cel.celldm) < INPUT.OO_distance/INPUT.unitconv)
    {
        for(int i=0; i<waters[O1].numH; i++)
        {
            pH = Cel.atoms[H].pos[waters[O1].indexH[i]];
            double angle = pO1.angle(pH,pO2,Cel.celldm);
            if( angle < INPUT.HOO_angle)
            {
                Hbonds[O1].donate_angle.push_back(angle);
                Hbonds[O1].donate_count++;
                Hbonds[O1].HB_count++;
                Hbonds[O1].index_donate.push_back(O2);
                Hbonds[O2].accept_count++;
                Hbonds[O2].HB_count++;
                Hbonds[O2].index_accept.push_back(O1);
                if(waters[O1].numH == 2)
                    avg_dwater++;
                else
                    avg_dion++;
                if(waters[O2].numH == 2)
                    avg_awater++;
                else
                    avg_aion++;
            }
        }
    }
}


