#include "reorganize.h"
#include "input.h"

void Reorganize::Routine()
{
    cout << "Reorganizing the System with type: " << INPUT.type << endl;
    switch (INPUT.type){
        case 0:
            cout << " Generate .xyz file normally!" << endl;
            break;
        case 1:
            cout << " Generate .xyz file with non-water highlighted!" << endl;
            break;
        case 2:
            cout << " Generate .xyz file with non-water only!" << endl;
            break;
        case 3:
            cout << " Generate the compatiable .pos and .wfc files" << endl;
            break;
        case 4:
            cout << " Generate the average value of NPT cell" << endl;
            break;
        default:
            cout << " Wrong INPUT type" << endl;
            break;
    }
    
    switch(INPUT.type){
        case 0:
        {
            ifstream ifs(INPUT.geo_file.c_str());
            file_assert(ifs,INPUT.geo_file);
            ofstream ofs1("new.xyz");
            file_assert(ofs1,"new.xyz");
            ofs1 << setprecision(10);
            
            cout << "Generating new.xyz file" << endl;
            
            for(int i=1; i<INPUT.ss_stop; i++)
            {
                Cellfile cel;
                
                if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
                {
                    if(cel.ReadGeometry(ifs,1)==0) break; // skip the one we dont need
                }
                else
                {
                    if(cel.ReadGeometry(ifs)==0) break; // if to the end of the file skip
                    
                    cel.organize_pos();
                    Reorganize::Print_xyz(ofs1,cel);
                    
                }
            }
        }
            break;
        case 1:
        {
            ifstream ifs(INPUT.geo_file.c_str());
            file_assert(ifs,INPUT.geo_file);
            ofstream ofs1("new_water.xyz");
            file_assert(ofs1,"new_water.xyz");
            ofs1 << setprecision(10);
            
            cout << "Generating new_water.xyz file" << endl;
            
            for(int i=1; i<INPUT.ss_stop; i++)
            {
                Cellfile cel;
                
                if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
                {
                    if(cel.ReadGeometry(ifs,1)==0) break; // skip the one we dont need
                }
                else
                {
                    if(cel.ReadGeometry(ifs)==0) break; // if to the end of the file skip
                    
                    cel.organize_pos();
                    Reorganize::Print_xyz_water(ofs1,cel);
                }
            }
        }
            break;
        case 2:
        {
            ifstream ifs(INPUT.geo_file.c_str());
            file_assert(ifs,INPUT.geo_file);
            ofstream ofs1("new_non_water.xyz");
            file_assert(ofs1,"new_non_water.xyz");
            ofs1 << setprecision(10);
            
            cout << "Generating new_non_water.xyz file" << endl;
            
            for(int i=1; i<INPUT.ss_stop; i++)
            {
                Cellfile cel;
                
                if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
                {
                    if(cel.ReadGeometry(ifs,1)==0) break; // skip the one we dont need
                }
                else
                {
                    if(cel.ReadGeometry(ifs)==0) break; // if to the end of the file skip
                    
                    cel.organize_pos();
                    Reorganize::Print_xyz_non_water(ofs1,cel);
                }
            }
        }
            break;
        case 3:
        {
            my_assert(INPUT.nband>0,"Please put in nband in INPUT");
            my_assert(INPUT.atom_num_tot>0,"Please put in atom information in INPUT");
            
            ifstream ifs_pos(INPUT.geo_file.c_str());
            file_assert(ifs_pos,INPUT.geo_file);
            ifstream ifs_wfc(INPUT.wan_file.c_str());
            file_assert(ifs_wfc,INPUT.wan_file);
            ofstream ofs_pos("new.pos");
            file_assert(ofs_pos,"new.pos");
            ofs_pos << setprecision(10);
            ofstream ofs_wfc("new.wfc");
            file_assert(ofs_wfc,"new.wfc");
            ofs_wfc << setprecision(10);
            
            int ss_pos,ss_wfc;
            double t_pos,t_wfc,t = 0;
            string line;
            
            ifs_pos >> ss_pos >> t_pos;
            ifs_pos.ignore(200,'\n');
            
            ifs_wfc >> ss_wfc >> t_wfc;
            ifs_wfc.ignore(200,'\n');
            
            while(ifs_pos.good() and ifs_wfc.good())
            {
                if(ss_pos == ss_wfc and abs(t_pos - t_wfc) < INPUT.eps)
                {
                    cout << "Writing  snapshot: " << ss_wfc << "  Time: " << t_wfc << "  dt= " << t_wfc - t << endl;
                    t = t_wfc;
                    ofs_wfc << ss_wfc << '\t' << t_wfc << endl;
                    for(int i=0; i<INPUT.nband; i++)
                    {
                        getline(ifs_wfc,line);
                        ofs_wfc << line << endl;
                    }
                    
                    ifs_wfc >> ss_wfc >> t_wfc;
                    ifs_wfc.ignore(200,'\n');
                    
                    ofs_pos << ss_pos << '\t' << t_pos << endl;
                    for(int i=0; i<INPUT.atom_num_tot; i++)
                    {
                        getline(ifs_pos,line);
                        ofs_pos << line << endl;
                    }
                    
                    ifs_pos >> ss_pos >> t_pos;
                    ifs_pos.ignore(200,'\n');
                }
                else if(ss_pos > ss_wfc and t_pos > t_wfc)
                {
                    cout << "Skipping snapshot: " << ss_wfc << " of .wfc"<< endl;
                    for(int i=0; i<INPUT.nband; i++)
                        ifs_wfc.ignore(200,'\n');
                    ifs_wfc >> ss_wfc >> t_wfc;
                    ifs_wfc.ignore(200,'\n');
                }
                else if(ss_pos < ss_wfc and t_pos < t_wfc)
                {
                    cout << "Skipping snapshot: " << ss_pos << " of .pos" << endl;
                    for(int i=0; i<INPUT.atom_num_tot; i++)
                        ifs_pos.ignore(200,'\n');
                    ifs_pos >> ss_pos >> t_pos;
                    ifs_pos.ignore(200,'\n');
                }
                else
                    my_assert(0,"The .pos and .wfc are not from the same calculation!");
            }
        }
            break;
        case 4:
        {
            Read_average_cel();
            
            ofstream ofs_cel("average.cel");
            file_assert(ofs_cel,"average.cel");
            
            ofs_cel << INPUT.Celldm << endl;
            cout << " Average celldm are\n" << INPUT.Celldm << endl;
        }
            break;
        default:
            break;
    }

}

void Reorganize::Print_xyz(ofstream &ofs, Cellfile &cel)
{
    assert(cel.allocate_atoms);
    
    ofs << cel.natom << endl;
    ofs << cel.snapshot << '\t' << cel.time << endl;
    
    for(int i=0; i<cel.ntype; i++)
    {
        for(int j=0; j<cel.atoms[i].na; j++)
        {
            ofs << cel.atoms[i].id << "  " << cel.atoms[i].pos[j].x*INPUT.unitconv << "  " << cel.atoms[i].pos[j].y*INPUT.unitconv << "  " << cel.atoms[i].pos[j].z*INPUT.unitconv << endl;
        }
    }
}

void Reorganize::Print_xyz_water(ofstream &ofs, Cellfile &cel)
{
    assert(cel.allocate_atoms);
    
    ofs << cel.natom << endl;
    ofs << cel.snapshot << '\t' << cel.time << endl;
    
    Waterfile WF;
    WF.Read_water(cel);
    
    for(int i=0; i<cel.ntype; i++)
    {
        for(int j=0; j<cel.atoms[i].na; j++)
        {
            if(i==WF.O && WF.waters[j].numH !=2 )
                ofs << "S  " << cel.atoms[i].pos[j].x*INPUT.unitconv << "  " << cel.atoms[i].pos[j].y*INPUT.unitconv << "  " << cel.atoms[i].pos[j].z*INPUT.unitconv << endl;
            else
                ofs << cel.atoms[i].id << "  " << cel.atoms[i].pos[j].x*INPUT.unitconv << "  " << cel.atoms[i].pos[j].y*INPUT.unitconv << "  " << cel.atoms[i].pos[j].z*INPUT.unitconv << endl;
        }
    }
}
    
void Reorganize::Print_xyz_non_water(ofstream &ofs, Cellfile &cel)
{
    assert(cel.allocate_atoms);
    
    Waterfile WF;
    WF.Read_water(cel);
    
    int natoms = 0;
    
    for(int i=0; i<WF.nion; i++)
    {
        natoms += WF.waters[WF.ions[i]].numH+1;
    }
    
    ofs << natoms << endl;
    ofs << cel.snapshot << '\t' << cel.time << endl;
    
    for(int i=0; i<WF.nion; i++)
    {
        ofs << "O  " << cel.atoms[WF.O].pos[WF.ions[i]]*INPUT.unitconv << endl;
        for(int j=0; j<WF.waters[WF.ions[i]].numH; j++)
            ofs << "H  " << cel.atoms[WF.H].pos[WF.waters[WF.ions[i]].indexH[j]]*INPUT.unitconv << endl;
    }
}

void Reorganize::Read_average_cel(int type)
{
    if(type == 1) cout << "Calculate average cell" << endl;
        
    my_assert(INPUT.ensemble == "npt", "Please make sure its npt system!");
    my_assert(INPUT.ifs_cel.good()," Can not open .cel file! ");
    
    INPUT.ifs_cel.seekg(INPUT.cel_top);
    int tmps,count=0;
    double tmpt;
    Vector3<double> tmpC;
    tmpC.set(0,0,0);
    
    for(int i=1; i<INPUT.ss_stop; i++)
    {
        if(INPUT.ifs_cel.eof()) break;
        
        INPUT.ifs_cel >> tmps >> tmpt;
        INPUT.ifs_cel.ignore(100,'\n');
        
        if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
        {   // skip the one we dont need
            if(type == 0)cout << "Skip Snapshot: " << tmps << " (cell)" << endl;
            Skip_lines(INPUT.ifs_cel,3);
        }
        else
        {
            if(type == 0)cout << "Read Snapshot: " << tmps << " (cell)" << endl;
            INPUT.ifs_cel >> INPUT.Celldm.x >> tmpt >> tmpt;
            INPUT.ifs_cel >> tmpt >> INPUT.Celldm.y >> tmpt;
            INPUT.ifs_cel >> tmpt >> tmpt >> INPUT.Celldm.z;
            count++;
            tmpC = tmpC + INPUT.Celldm;
        }
    }
    tmpC /= count;
    INPUT.Celldm = tmpC;
    INPUT.ss_n = count;
}
