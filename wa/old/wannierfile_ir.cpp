#include "wannierfile.h"
#include "input.h"
#include "reorganize.h"

void Wannierfile::Routine_Ir()
{
    cout << "Using Wannier Center Calculating IR spectra using type: " << INPUT.type << endl;
    switch (INPUT.type) {
        case 0:
            cout << "Start from Reading .wfc .pos (.vel) files" << endl;
            break;
        case 1:
            cout << "Generate dipole.txt and vdipole.txt only" << endl;
            break;
        case 2:
            cout << "Start from Reading vdipole.txt file" << endl;
            break;
        case 3:
            cout << "Calculate Intramolecular contribution from vdipole.txt file" << endl;
            break;
        case 31:
            cout << "Calculate Spatial Intermolecular contribution (pos_start0 to pos_start1) from vdipole.txt file" << endl;
            break;
        case 311:
            cout << "Calculate Spatial Intermolecular contribution (pos_start0 to pos_start1) from dipole.txt and vdipole.txt file for angle>0" << endl;
            break;
        case 32:
            cout << "Calculate Spatial Intermolecular contribution (pos_start0 to pos_start1) from vdipole.txt file to several segments at the same time" << endl;
            break;
        case 4:
            cout << "Calculate electronic and ionic contributions" << endl;
            break;
        case 41:
            cout << "Calculate lone-pair and bond-pair electronic contributions" << endl;
            break;
        case 5:
            cout << "Generate electronic and ionic contributions to edipole.txt and idipole.txt only" << endl;
            break;
        case 51:
            cout << "Generate lone-pair and bond-pair electronic contributions to eldipole.txt and ebdipole.txt only" << endl;
            break;
        case 6:
            cout << "Calculate electronic and ionic contributions from edipole.txt and idipole.txt" << endl;
            break;
        case 7:
            cout << "Calculate elctronic ionic coupling term from ion,electron_vdipole.txt" << endl;
            break;
        case 71:
            cout << "Calculate lone-pair and bond-pair electronic coupling term from electron_lone/bond_vdipole.txt" << endl;
            break;
        default:
            cout << "Wrong type of IR" << endl;
            break;
    }
    
    cout << setprecision(10);
    
    int tjump=0;
    Vector3<double> *vdipole_file; // vdipole[snapshot j] (Debye/a.u.time)
    if(INPUT.type == 0 or INPUT.type == 2 or INPUT.type == 3 or INPUT.type == 4 or INPUT.type == 41 or INPUT.type == 6)
        vdipole_file = new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
    
    int O;
    for(int i; i< INPUT.ntype; i++)
        if(INPUT.atom_name[i] == "O")
            O = i;
    
    switch (INPUT.type) {
        case 0:
        {
            ofstream ofs("waters_dipole.txt");
            ofstream ofs_v("waters_vdipole.txt");
            Wannierfile::Routine_Ir_1(ofs,ofs_v,vdipole_file);
            //generate tcf to list
            double *tcf;
            tcf = new double [INPUT.delta];
            for(int i=0; i<INPUT.delta; i++)
                tcf[i] = 0;
            TCF2list(tcf,vdipole_file,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha);
            //write TCF to output data
            ofstream ofs_tcf("TCF.txt");
            file_assert(ofs_tcf,"TCF.txt");
            ofs_tcf << setprecision(10);
            for(int i=0;i<INPUT.delta;i++)
                ofs_tcf << INPUT.dt*i*(1+tjump) << "  " << tcf[i] << endl;
            //FT(tcf,INPUT.dt,INPUT.delta-1,INPUT.delta);
            FT(tcf,INPUT.dt*(1+tjump),INPUT.delta-1);
        }
            break;
        case 1:
        {
            ofstream ofs("waters_dipole.txt");
            ofstream ofs_v("waters_vdipole.txt");
            Wannierfile::Routine_Ir_1(ofs,ofs_v);
        }
            break;
        case 2:
        {
            assert(INPUT.dt != -1.0);
            
            if(INPUT.in_file == "none")
                INPUT.in_file = "waters_vdipole.txt";
            ifstream ifs(INPUT.in_file.c_str());
            file_assert(ifs,INPUT.in_file);
            INPUT.ss_n = Wannierfile::Read_Vdipolefile(ifs,vdipole_file,INPUT.atom_num[O]);
            //generate tcf to list
            double *tcf;
            tcf = new double [INPUT.delta];
            for(int i=0; i<INPUT.delta; i++)
                tcf[i] = 0;
            TCF2list(tcf,vdipole_file,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha);
            //write TCF to output data
            ofstream ofs_tcf("TCF.txt");
            file_assert(ofs_tcf,"TCF.txt");
            ofs_tcf << setprecision(10);
            for(int i=0;i<INPUT.delta;i++)
                ofs_tcf << INPUT.dt*i*(1+tjump) << "  " << tcf[i] << endl;
            //FT(tcf,INPUT.dt,INPUT.delta-1,INPUT.delta);
            FT(tcf,INPUT.dt*(1+tjump),INPUT.delta-1);
        }
            break;
        case 3:
        {
            assert(INPUT.dt != -1.0);
            if(INPUT.in_file == "none")
                INPUT.in_file = "waters_vdipole.txt";
            
            int O;
            for(int i; i< INPUT.ntype; i++)
                if(INPUT.atom_name[i] == "O")
                    O = i;
            //generate a list to store the result of TCF calculation
            double *TCF_result;
            int count=0;
            TCF_result = new double [INPUT.delta];
            for(int i=0; i<INPUT.delta; i++)
                TCF_result[i] = 0;
            //read vdipoles in [water i][snapshot j]
            Vector3<double> **vdipole_test;
            vdipole_test = new Vector3<double> * [INPUT.atom_num[O]];
            for(int i=0; i<INPUT.atom_num[O];i++)
                vdipole_test[i] =  new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            
            ifstream ifs(INPUT.in_file.c_str());
            file_assert(ifs,INPUT.in_file);
            INPUT.ss_n = Wannierfile::Read_Vdipolefile(ifs,vdipole_test,INPUT.atom_num[O]);
            //calculate only intracontribution
            for(int i=0; i<INPUT.atom_num[O];i++)
            {
                int tmp = TCF2list(TCF_result,vdipole_test[i],INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha,i==INPUT.atom_num[O]-1);
                if(tmp > count)
                    count = tmp;
            }
            //print out the intracontribution
            ofstream ofs_intra("TCF_intra.txt");
            ofs_intra << setprecision(10);
            for(int i=0; i<count; i++)
                ofs_intra << INPUT.dt*i*(1+tjump) << " " << TCF_result[i] << endl;
            //generate IR spectra
            FT(TCF_result,INPUT.dt*(1+tjump),count-1,"FT_intra.txt");
            //generate FT_inter.txt if FT.txt and FT_intra.txt exist
            cout << "Generate FT_inter.txt if possible " << endl;
            ifstream ifs_FT("FT.txt");
            ifstream ifs_FT_intra("FT_intra.txt");
            ofstream ofs_FT_inter("FT_inter.txt");
            double time1,time2,FT1,FT2;
            while(ifs_FT.good())
            {
                ifs_FT >> time1;
                ifs_FT_intra >> time2;
                if(ifs_FT_intra.good())
                {
                    if(time1 != time2)
                    {
                        cout << "FT.txt and FT_intra.txt are not Compatable, Stop" << endl;
                        exit(1);
                    }
                    Read_value(ifs_FT,FT1);
                    Read_value(ifs_FT_intra,FT2);
                    ofs_FT_inter << time1 << " " << FT1-FT2 << endl;
                }
            }
        }
            break;
        case 31:
        {
            my_assert(INPUT.pos_start[0] < INPUT.pos_start[1], "Please give pos_start[1][2] (r1,r2) in INPUT!");
            my_assert(INPUT.pos_start_n != 0, "Please give pos_start_n in INPUT!");
            my_assert(INPUT.dt != -1.0, "Please give dt in INPUT!");
            if(INPUT.in_file == "none")
                INPUT.in_file = "waters_vdipole.txt";
            
            int O;
            for(int i; i< INPUT.ntype; i++)
                if(INPUT.atom_name[i] == "O")
                    O = i;
            
            // generate dipole_pos[water i][snapshot j] (SI)
            // generate celldm[snapshot j] (SI)
            // generate vdipole[water i][snapshot j] (Debye/a.u.time)
            cout << "Generate dipole_pos and celldm and vdipole" << endl;
            INPUT.ifs_cel.seekg(INPUT.cel_top);
            ifstream ifs_geo(INPUT.geo_file.c_str());
            file_assert(ifs_geo,INPUT.geo_file);
            ifstream ifs(INPUT.in_file.c_str());
            file_assert(ifs,INPUT.in_file);
            
            cout << "Reading vdipolefile " << endl;
            Vector3<double> **vdipole_test;
            vdipole_test = new Vector3<double> * [INPUT.atom_num[O]];
            for(int i=0; i<INPUT.atom_num[O];i++)
                vdipole_test[i] =  new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            INPUT.ss_n = Wannierfile::Read_Vdipolefile(ifs,vdipole_test,INPUT.atom_num[O]);
            
            Vector3<double> **dipole_pos;
            dipole_pos = new Vector3<double> * [INPUT.atom_num[O]];
            for(int i=0; i<INPUT.atom_num[O];i++)
                dipole_pos[i] =  new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            
            Vector3<double> *celldm;
            celldm = new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            
            int count = 0;
            
            cout << "Reading celldm and dipole_pos" << endl;
            for(int i=1; i<INPUT.ss_stop; i++)
            {
                Cellfile cel;
                
                if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
                {
                    cel.ReadGeometry(ifs_geo,1); // skip the one we dont need
                }
                else
                {
                    if(cel.ReadGeometry(ifs_geo)==0) break; // if to the end of the file skip
                    for(int j=0; j<cel.atoms[O].na; j++)
                        dipole_pos[j][count] = cel.atoms[O].pos[j]*INPUT.unitconv;
                    celldm[count] = cel.celldm*INPUT.unitconv;
                    count ++;
                }
            }
            
            count -= 2;
            //see if dipole_pos and vdipole compatable
            my_assert(INPUT.ss_n == count,"dipole_pos and vdipole are not compatable!");
            
            //generate a list to store the result of TCF calculation
            double *TCF_result;
            TCF_result = new double [INPUT.delta];
            
            //calculate TCF from pos_start[0] to [1] with pos_start_n segments
            double start = INPUT.pos_start[0], end = INPUT.pos_start[1], step = (end - start)/INPUT.pos_start_n;
            end = start + step;
            for(int k=0; k<INPUT.pos_start_n; k++)
            {
                for(int i=0; i<INPUT.delta; i++)
                    TCF_result[i] = 0;
                count = 0;
                for(int i=0; i<INPUT.atom_num[O];i++)
                {
                    for(int j=0; j< INPUT.atom_num[O]; j++)
                    {
                        if(i!=j)
                        {
                            int tmp = TCF2list3(TCF_result,vdipole_test[i],vdipole_test[j],dipole_pos[i],dipole_pos[j],celldm,start,end,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha,count==0);
                            if(tmp > count)
                                count = tmp;
                            cout << "Correlating vdipole " << i << " and vdipole " << j << " from r1= " << start << " r2= " << end << endl;
                        }
                    }
                }
                //print out the intracontribution
                ostringstream ss1,ss2;
                ss1 << start;
                ss2 << end;
                string tmp1(ss1.str());
                string tmp2(ss2.str());
                ofstream ofs_inter_spatial(("TCF_inter_spatial_" + tmp1 + "_" + tmp2 + ".txt").c_str());
                ofs_inter_spatial << setprecision(10);
                for(int i=0; i<count; i++)
                    ofs_inter_spatial << INPUT.dt*i*(1+tjump) << " " << TCF_result[i] << endl;
                //generate IR spectra
                FT(TCF_result,INPUT.dt*(1+tjump),count-1,"FT_inter_spatial_" + tmp1 + "_" + tmp2 + ".txt");
                start += step;
                end += step;
            }
        }
            break;
        case 311:
        {
            my_assert(INPUT.pos_start[0] < INPUT.pos_start[1], "Please give pos_start[1][2] (r1,r2) in INPUT!");
            my_assert(INPUT.pos_start_n != 0, "Please give pos_start_n in INPUT!");
            my_assert(INPUT.dt != -1.0, "Please give dt in INPUT!");
            if(INPUT.in_file == "none")
                INPUT.in_file = "./";
            
            int O;
            for(int i; i< INPUT.ntype; i++)
                if(INPUT.atom_name[i] == "O")
                    O = i;
            
            // generate dipole_pos[water i][snapshot j] (SI)
            // generate celldm[snapshot j] (SI)
            // generate vdipole[water i][snapshot j] (Debye/a.u.time)
            cout << "Generate dipole_pos and celldm and vdipole" << endl;
            INPUT.ifs_cel.seekg(INPUT.cel_top);
            ifstream ifs_geo(INPUT.geo_file.c_str());
            file_assert(ifs_geo,INPUT.geo_file);
            ifstream ifs((INPUT.in_file+"waters_vdipole.txt").c_str());
            file_assert(ifs,INPUT.in_file+"waters_vdipole.txt");
            ifstream ifs2((INPUT.in_file+"waters_dipole.txt").c_str());
            file_assert(ifs2,INPUT.in_file+"waters_dipole.txt");
            
            cout << "Reading vdipolefile " << endl;
            Vector3<double> **vdipole_test;
            vdipole_test = new Vector3<double> * [INPUT.atom_num[O]];
            for(int i=0; i<INPUT.atom_num[O];i++)
                vdipole_test[i] =  new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            INPUT.ss_n = Wannierfile::Read_Vdipolefile(ifs,vdipole_test,INPUT.atom_num[O]);
            
            Vector3<double> **dipole_test;
            dipole_test = new Vector3<double> * [INPUT.atom_num[O]];
            for(int i=0; i<INPUT.atom_num[O];i++)
                dipole_test[i] =  new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            Wannierfile::Read_Dipolefile(ifs2,dipole_test,INPUT.atom_num[O]);
            
            Vector3<double> **dipole_pos;
            dipole_pos = new Vector3<double> * [INPUT.atom_num[O]];
            for(int i=0; i<INPUT.atom_num[O];i++)
                dipole_pos[i] =  new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            
            Vector3<double> *celldm;
            celldm = new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            
            int count = 0;
            
            cout << "Reading celldm and dipole_pos" << endl;
            for(int i=1; i<INPUT.ss_stop; i++)
            {
                Cellfile cel;
                
                if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
                {
                    cel.ReadGeometry(ifs_geo,1); // skip the one we dont need
                }
                else
                {
                    if(cel.ReadGeometry(ifs_geo)==0) break; // if to the end of the file skip
                    for(int j=0; j<cel.atoms[O].na; j++)
                        dipole_pos[j][count] = cel.atoms[O].pos[j]*INPUT.unitconv;
                    celldm[count] = cel.celldm*INPUT.unitconv;
                    count ++;
                }
            }
            
            count -= 2;
            //see if dipole_pos and vdipole compatable
            my_assert(INPUT.ss_n == count,"dipole_pos and vdipole are not compatable!");
            
            //generate a list to store the result of TCF calculation
            double *TCF_result;
            TCF_result = new double [INPUT.delta];
            
            //calculate TCF from pos_start[0] to [1] with pos_start_n segments
            double start = INPUT.pos_start[0], end = INPUT.pos_start[1], step = (end - start)/INPUT.pos_start_n;
            end = start + step;
            for(int k=0; k<INPUT.pos_start_n; k++)
            {
                for(int i=0; i<INPUT.delta; i++)
                    TCF_result[i] = 0;
                count = 0;
                for(int i=0; i<INPUT.atom_num[O];i++)
                {
                    for(int j=0; j< INPUT.atom_num[O]; j++)
                    {
                        if(i!=j)
                        {
                            int tmp = TCF2list5(TCF_result,vdipole_test[i],vdipole_test[j],dipole_test[i],dipole_test[j],dipole_pos[i],dipole_pos[j],celldm,start,end,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha,count==0);
                            if(tmp > count)
                                count = tmp;
                            cout << "Correlating vdipole " << i << " and vdipole " << j << " from r1= " << start << " r2= " << end << endl;
                        }
                    }
                }
                //print out the intracontribution
                ostringstream ss1,ss2;
                ss1 << start;
                ss2 << end;
                string tmp1(ss1.str());
                string tmp2(ss2.str());
                ofstream ofs_inter_spatial(("TCF_inter_spatial_" + tmp1 + "_" + tmp2 + ".txt").c_str());
                ofs_inter_spatial << setprecision(10);
                for(int i=0; i<count; i++)
                    ofs_inter_spatial << INPUT.dt*i*(1+tjump) << " " << TCF_result[i] << endl;
                //generate IR spectra
                FT(TCF_result,INPUT.dt*(1+tjump),count-1,"FT_inter_spatial_" + tmp1 + "_" + tmp2 + ".txt");
                start += step;
                end += step;
            }
        }
            break;
        case 32:
        {
            my_assert(INPUT.pos_start[0] < INPUT.pos_start[1], "Please give pos_start[1][2] (r1,r2) in INPUT!");
            my_assert(INPUT.pos_start_n != 0, "Please give pos_start_n in INPUT!");
            my_assert(INPUT.dt != -1.0, "Please give dt in INPUT!");
            if(INPUT.in_file == "none")
                INPUT.in_file = "waters_vdipole.txt";
            
            int O;
            for(int i; i< INPUT.ntype; i++)
                if(INPUT.atom_name[i] == "O")
                    O = i;
            
            // generate dipole_pos[water i][snapshot j] (SI)
            // generate celldm[snapshot j] (SI)
            // generate vdipole[water i][snapshot j] (Debye/a.u.time)
            cout << "Generate dipole_pos and celldm and vdipole" << endl;
            INPUT.ifs_cel.seekg(INPUT.cel_top);
            ifstream ifs_geo(INPUT.geo_file.c_str());
            file_assert(ifs_geo,INPUT.geo_file);
            ifstream ifs(INPUT.in_file.c_str());
            file_assert(ifs,INPUT.in_file);
            
            cout << "Reading vdipolefile " << endl;
            Vector3<double> **vdipole_test;
            vdipole_test = new Vector3<double> * [INPUT.atom_num[O]];
            for(int i=0; i<INPUT.atom_num[O];i++)
                vdipole_test[i] =  new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            INPUT.ss_n = Wannierfile::Read_Vdipolefile(ifs,vdipole_test,INPUT.atom_num[O]);
            
            Vector3<double> **dipole_pos;
            dipole_pos = new Vector3<double> * [INPUT.atom_num[O]];
            for(int i=0; i<INPUT.atom_num[O];i++)
                dipole_pos[i] =  new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            
            Vector3<double> *celldm;
            celldm = new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            
            int count = 0;
            
            cout << "Reading celldm and dipole_pos" << endl;
            for(int i=1; i<INPUT.ss_stop; i++)
            {
                Cellfile cel;
                
                if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
                {
                    cel.ReadGeometry(ifs_geo,1); // skip the one we dont need
                }
                else
                {
                    if(cel.ReadGeometry(ifs_geo)==0) break; // if to the end of the file skip
                    for(int j=0; j<cel.atoms[O].na; j++)
                        dipole_pos[j][count] = cel.atoms[O].pos[j]*INPUT.unitconv;
                    celldm[count] = cel.celldm*INPUT.unitconv;
                    count ++;
                }
            }
            
            count -= 2;
            //see if dipole_pos and vdipole compatable
            my_assert(INPUT.ss_n == count,"dipole_pos and vdipole are not compatable!");
            
            //generate a list to store the result of TCF calculation
            double **TCF_result;
            TCF_result = new double * [INPUT.pos_start_n];
            for(int i=0; i<INPUT.pos_start_n; i++)
                TCF_result[i] = new double [INPUT.delta];
            
            //calculate TCF from pos_start[0] to [1] with pos_start_n segments
            double start = INPUT.pos_start[0], end = INPUT.pos_start[1], step = (end - start)/INPUT.pos_start_n;
            for(int i=0; i<INPUT.pos_start_n; i++)
                for(int j=0; j<INPUT.delta; j++)
                    TCF_result[i][j] = 0;
            count = 0;
            for(int i=0; i<INPUT.atom_num[O];i++)
            {
                for(int j=0; j< INPUT.atom_num[O]; j++)
                {
                    if(i!=j)
                    {
                        int tmp = TCF2list4(TCF_result,vdipole_test[i],vdipole_test[j],dipole_pos[i],dipole_pos[j],celldm,start,end,INPUT.pos_start_n,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha,count==0);
                        if(tmp > count)
                            count = tmp;
                        cout << "Correlating vdipole " << i << " and vdipole " << j << endl;
                    }
                }
            }
            //print out the intracontribution
            ostringstream ss1,ss2,ss3;
            ss1 << start;
            ss2 << end;
            ss3 << step;
            string tmp1(ss1.str());
            string tmp2(ss2.str());
            string tmp3(ss3.str());
            ofstream ofs_inter_spatial(("TCF_inter_spatial_" + tmp1 + "_" + tmp2 + "_step_" + tmp3 + ".txt").c_str());
            ofs_inter_spatial << setprecision(10);
            for(int i=0; i<count; i++)
            {
                ofs_inter_spatial << INPUT.dt*i*(1+tjump);
                for(int j=0; j<INPUT.pos_start_n ; j++)
                    ofs_inter_spatial << " " << TCF_result[j][i];
                ofs_inter_spatial << endl;
            }
            
            //generate IR spectra
            for(int i=0; i<INPUT.pos_start_n; i++)
            {
                ostringstream ss;
                ss << i;
                string tmp(ss.str());
                FT(TCF_result[i],INPUT.dt*(1+tjump),count-1,"FT_inter_spatial_" + tmp1 + "_" + tmp2 + "_step_" + tmp3 + "_" + tmp +".txt");
            }
        }
            break;
        case 4:
        {
            {
                ofstream ofs("ion_dipole.txt");
                ofstream ofs_v("ion_vdipole.txt");
                Wannierfile::Routine_Ir_1(ofs,ofs_v,vdipole_file,1);
                //generate tcf to list
                double *tcf;
                tcf = new double [INPUT.delta];
                for(int i=0; i<INPUT.delta; i++)
                    tcf[i] = 0; TCF2list(tcf,vdipole_file,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha);
                //write TCF to output data
                ofstream ofs_tcf("TCF_ion.txt");
                file_assert(ofs_tcf,"TCF_ion.txt");
                ofs_tcf << setprecision(10);
                for(int i=0;i<INPUT.delta;i++)
                    ofs_tcf << INPUT.dt*i*(1+tjump) << "  " << tcf[i] << endl;
                //FT(tcf,INPUT.dt,INPUT.delta-1,INPUT.delta);
                FT(tcf,INPUT.dt*(1+tjump),INPUT.delta-1,"FT_ion.txt");
            }
            {
                ofstream ofs("electron_dipole.txt");
                ofstream ofs_v("electron_vdipole.txt");
                Wannierfile::Routine_Ir_1(ofs,ofs_v,vdipole_file,2);
                //generate tcf to list
                double *tcf;
                tcf = new double [INPUT.delta];
                for(int i=0; i<INPUT.delta; i++)
                    tcf[i] = 0;
                TCF2list(tcf,vdipole_file,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha);
                //write TCF to output data
                ofstream ofs_tcf("TCF_electron.txt");
                file_assert(ofs_tcf,"TCF_electron.txt");
                ofs_tcf << setprecision(10);
                for(int i=0;i<INPUT.delta;i++)
                    ofs_tcf << INPUT.dt*i*(1+tjump) << "  " << tcf[i] << endl;
                //FT(tcf,INPUT.dt,INPUT.delta-1,INPUT.delta);
                FT(tcf,INPUT.dt*(1+tjump),INPUT.delta-1,"FT_electron.txt");
            }
        }
            break;
        case 41:
        {
            {
                ofstream ofs("electron_lone_dipole.txt");
                ofstream ofs_v("electron_lone_vdipole.txt");
                Wannierfile::Routine_Ir_1(ofs,ofs_v,vdipole_file,3);
                //generate tcf to list
                double *tcf;
                tcf = new double [INPUT.delta];
                for(int i=0; i<INPUT.delta; i++)
                    tcf[i] = 0;
                TCF2list(tcf,vdipole_file,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha);
                //write TCF to output data
                ofstream ofs_tcf("TCF_electron_lone.txt");
                file_assert(ofs_tcf,"TCF_electron_lone.txt");
                ofs_tcf << setprecision(10);
                for(int i=0;i<INPUT.delta;i++)
                    ofs_tcf << INPUT.dt*i*(1+tjump) << "  " << tcf[i] << endl;
                //FT(tcf,INPUT.dt,INPUT.delta-1,INPUT.delta);
                FT(tcf,INPUT.dt*(1+tjump),INPUT.delta-1,"FT_electron_lone.txt");
            }
            {
                ofstream ofs("electron_bond_dipole.txt");
                ofstream ofs_v("electron_bond_vdipole.txt");
                Wannierfile::Routine_Ir_1(ofs,ofs_v,vdipole_file,4);
                //generate tcf to list
                double *tcf;
                tcf = new double [INPUT.delta];
                for(int i=0; i<INPUT.delta; i++)
                    tcf[i] = 0;
                TCF2list(tcf,vdipole_file,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha);
                //write TCF to output data
                ofstream ofs_tcf("TCF_electron_bond.txt");
                file_assert(ofs_tcf,"TCF_electron_bond.txt");
                ofs_tcf << setprecision(10);
                for(int i=0;i<INPUT.delta;i++)
                    ofs_tcf << INPUT.dt*i*(1+tjump) << "  " << tcf[i] << endl;
                //FT(tcf,INPUT.dt,INPUT.delta-1,INPUT.delta);
                FT(tcf,INPUT.dt*(1+tjump),INPUT.delta-1,"FT_electron_bond.txt");
            }
        }
            break;
        case 5:
        {
            ofstream ofs1("ion_dipole.txt");
            ofstream ofs1_v("ion_vdipole.txt");
            ofstream ofs2("electron_dipole.txt");
            ofstream ofs2_v("electron_vdipole.txt");
            
            INPUT.ifs_cel.seekg(INPUT.cel_top);
            ifstream ifs(INPUT.geo_file.c_str());
            ifstream ifs_wan(INPUT.wan_file.c_str());
            file_assert(ifs,INPUT.geo_file);
            file_assert(ifs_wan,INPUT.wan_file);
            ofs1 << setprecision(10);
            ofs1_v << setprecision(10);
            ofs2 << setprecision(10);
            ofs2_v << setprecision(10);
            
            Wannierfile p2_WF1;
            Wannierfile p_WF1;
            Wannierfile p2_WF2;
            Wannierfile p_WF2;
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
                    if(cel.ReadWC(ifs_wan)==0)break;
                    cel.organize_pos();
                    cel.organize_wan();
                    //read and write water_dipole
                    Wannierfile WF;
                    WF.Read_Wannier(cel);
                    WF.Read_Dipole(cel,1);
                    WF.Print_Dipole(ofs1);
                    //recort dt
                    if(INPUT.ss_n++ == 1)
                        INPUT.dt = WF.time - p_WF1.time;
                    //read and write water_vdipole and record tot_vdipole of each snapshot
                    WF.Read_Vdipole(p_WF1,p2_WF1,ofs1_v);
                    
                    WF.Read_Dipole(cel,2);
                    WF.Print_Dipole(ofs2);
                    WF.Read_Vdipole(p_WF2,p2_WF2,ofs2_v);
                }
            }
            cout << " dt: " << INPUT.dt << " total snapshot: "<< INPUT.ss_n << endl;
        }
            break;
        case 51:
        {
            ofstream ofs1("electron_lone_dipole.txt");
            ofstream ofs1_v("electron_lone_vdipole.txt");
            ofstream ofs2("electron_bond_dipole.txt");
            ofstream ofs2_v("electron_bond_vdipole.txt");
            
            INPUT.ifs_cel.seekg(INPUT.cel_top);
            ifstream ifs(INPUT.geo_file.c_str());
            ifstream ifs_wan(INPUT.wan_file.c_str());
            file_assert(ifs,INPUT.geo_file);
            file_assert(ifs_wan,INPUT.wan_file);
            ofs1 << setprecision(10);
            ofs1_v << setprecision(10);
            ofs2 << setprecision(10);
            ofs2_v << setprecision(10);
            
            Wannierfile p2_WF1;
            Wannierfile p_WF1;
            Wannierfile p2_WF2;
            Wannierfile p_WF2;
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
                    if(cel.ReadWC(ifs_wan)==0)break;
                    cel.organize_pos();
                    cel.organize_wan();
                    //read and write water_dipole
                    Wannierfile WF;
                    WF.Read_Wannier(cel);
                    WF.Read_Dipole(cel,3);
                    WF.Print_Dipole(ofs1);
                    //recort dt
                    if(INPUT.ss_n++ == 1)
                        INPUT.dt = WF.time - p_WF1.time;
                    //read and write water_vdipole and record tot_vdipole of each snapshot
                    WF.Read_Vdipole(p_WF1,p2_WF1,ofs1_v);
                    
                    WF.Read_Dipole(cel,4);
                    WF.Print_Dipole(ofs2);
                    WF.Read_Vdipole(p_WF2,p2_WF2,ofs2_v);
                }
            }
            cout << " dt: " << INPUT.dt << " total snapshot: "<< INPUT.ss_n << endl;
        }
            break;
        case 6:
        {
            assert(INPUT.dt != -1.0);
            {
                INPUT.in_file = "ion_vdipole.txt";
                ifstream ifs(INPUT.in_file.c_str());
                file_assert(ifs,INPUT.in_file);
                INPUT.ss_n = Wannierfile::Read_Vdipolefile(ifs,vdipole_file,INPUT.atom_num[O]);
                //generate tcf to list
                double *tcf;
                tcf = new double [INPUT.delta];
                for(int i=0; i<INPUT.delta; i++)
                    tcf[i] = 0;
                TCF2list(tcf,vdipole_file,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha);
                //write TCF to output data
                ofstream ofs_tcf("TCF_ion.txt");
                file_assert(ofs_tcf,"TCF_ion.txt");
                ofs_tcf << setprecision(10);
                for(int i=0;i<INPUT.delta;i++)
                    ofs_tcf << INPUT.dt*i*(1+tjump) << "  " << tcf[i] << endl;
                //FT(tcf,INPUT.dt,INPUT.delta-1,INPUT.delta);
                FT(tcf,INPUT.dt*(1+tjump),INPUT.delta-1,"FT_ion.txt");
            }
            {
                INPUT.in_file = "electron_vdipole.txt";
                ifstream ifs(INPUT.in_file.c_str());
                file_assert(ifs,INPUT.in_file);
                INPUT.ss_n = Wannierfile::Read_Vdipolefile(ifs,vdipole_file,INPUT.atom_num[O]);
                //generate tcf to list
                double *tcf;
                tcf = new double [INPUT.delta];
                for(int i=0; i<INPUT.delta; i++)
                    tcf[i] = 0;
                TCF2list(tcf,vdipole_file,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha);
                //write TCF to output data
                ofstream ofs_tcf("TCF_electron.txt");
                file_assert(ofs_tcf,"TCF_electron.txt");
                ofs_tcf << setprecision(10);
                for(int i=0;i<INPUT.delta;i++)
                    ofs_tcf << INPUT.dt*i*(1+tjump) << "  " << tcf[i] << endl;
                //FT(tcf,INPUT.dt,INPUT.delta-1,INPUT.delta);
                FT(tcf,INPUT.dt*(1+tjump),INPUT.delta-1,"FT_water.txt");
            }
        }
            break;
        case 7:
        {
            assert(INPUT.dt != -1);
            ifstream ifs_i("ion_vdipole.txt");
            file_assert(ifs_i,"Can not find ion_vdipole\n Please generate it using method 5",1);
            ifstream ifs_e("electron_vdipole.txt");
            file_assert(ifs_e,"Can not file electron_vdipole\n Please generate it using method 5",1);
            
            Vector3<double> *vdipole_file_i,*vdipole_file_e;
            vdipole_file_i = new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            vdipole_file_e = new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            
            INPUT.ss_n = Wannierfile::Read_Vdipolefile(ifs_i,vdipole_file_i,INPUT.atom_num[O]);
            my_assert(INPUT.ss_n == Wannierfile::Read_Vdipolefile(ifs_e,vdipole_file_e,INPUT.atom_num[O]),"ion_vdipole.txt and electron_vdipole.txt are not Compatable, Stop");
            
            //for(int i=0; i<INPUT.ss_n; i++)
            // vdipole_file_i[i].print2screen();
            
            //generate tcf to list
            double *tcf;
            tcf = new double [INPUT.delta];
            for(int i=0; i<INPUT.delta; i++)
                tcf[i] = 0;
            TCF2list2(tcf,vdipole_file_e,vdipole_file_i,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha,0);
            TCF2list2(tcf,vdipole_file_i,vdipole_file_e,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha);
            //write TCF to output data
            ofstream ofs_tcf("TCF_ie_couple.txt");
            file_assert(ofs_tcf,"TCF_ie_couple.txt");
            ofs_tcf << setprecision(10);
            for(int i=0;i<INPUT.delta;i++)
                ofs_tcf << INPUT.dt*i*(1+tjump) << "  " << tcf[i] << endl;
            FT(tcf,INPUT.dt*(1+tjump),INPUT.delta-1,"FT_ie_couple.txt");
        }
            break;
        case 71:
        {
            assert(INPUT.dt != -1);
            ifstream ifs_i("electron_lone_vdipole.txt");
            file_assert(ifs_i,"Can not find electron_lone_vdipole\n Please generate it using method 5",1);
            ifstream ifs_e("electron_bond_vdipole.txt");
            file_assert(ifs_e,"Can not file electron_bond_vdipole\n Please generate it using method 5",1);
            
            Vector3<double> *vdipole_file_i,*vdipole_file_e;
            vdipole_file_i = new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            vdipole_file_e = new Vector3<double> [(INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1];
            
            INPUT.ss_n = Wannierfile::Read_Vdipolefile(ifs_i,vdipole_file_i,INPUT.atom_num[O]);
            my_assert(INPUT.ss_n == Wannierfile::Read_Vdipolefile(ifs_e,vdipole_file_e,INPUT.atom_num[O]),"electron_lone_vdipole.txt and electron_bond_vdipole.txt are not Compatable, Stop");
            
            //for(int i=0; i<INPUT.ss_n; i++)
            // vdipole_file_i[i].print2screen();
            
            //generate tcf to list
            double *tcf;
            tcf = new double [INPUT.delta];
            for(int i=0; i<INPUT.delta; i++)
                tcf[i] = 0;
            TCF2list2(tcf,vdipole_file_e,vdipole_file_i,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha,0);
            TCF2list2(tcf,vdipole_file_i,vdipole_file_e,INPUT.dt,INPUT.ss_n,INPUT.delta,tjump,INPUT.cutoff,INPUT.alpha);
            //write TCF to output data
            ofstream ofs_tcf("TCF_lb_electron_couple.txt");
            file_assert(ofs_tcf,"TCF_lb_electron_couple.txt");
            ofs_tcf << setprecision(10);
            for(int i=0;i<INPUT.delta;i++)
                ofs_tcf << INPUT.dt*i*(1+tjump) << "  " << tcf[i] << endl;
            FT(tcf,INPUT.dt*(1+tjump),INPUT.delta-1,"FT_lb_electron_couple.txt");
        }
            /*case 8:
             {
             //read vdipole file
             if(INPUT.in_file == "none")
             INPUT.in_file = "water_vdipole.txt";
             ifstream ifs_vdipole(INPUT.in_file.c_str());
             file_assert(ifs_vdipole,INPUT.in_file);
             
             ifstream ifs(INPUT.geo_file.c_str());
             file_assert(ifs,INPUT.geo_file);
             
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
             //read water
             Waterfile WF;
             WF.Read_water();
             
             if(INPUT.ss_n != 0)
             {
             int x = INPUT.ss_n - 1;
             vdipole_file[x].set(0,0,0)
             for(int j=0; j<WF.nwater; j++)
             {
             if(WF.numH == 2)
             {
             Vector3<double> tmp;
             ifs_vdipole >> tmp.x >> tmp.y >> tmp.z;
             vdipole_file[x] = vdipole_file[x] + tmp*();
             }
             }
             }
             INPUT.ss_n++ == 1;
             //read and write water_vdipole and record tot_vdipole of each snapshot
             }
             }
             }*/
        default:
            break;
    }
}

void Wannierfile::Routine_Ir_1(ofstream &ofs, ofstream &ofs_v, int read_type)
{
    INPUT.ifs_cel.seekg(INPUT.cel_top);
    ifstream ifs(INPUT.geo_file.c_str());
    ifstream ifs_wan(INPUT.wan_file.c_str());
    file_assert(ifs,INPUT.geo_file);
    file_assert(ifs_wan,INPUT.wan_file);
    ofs << setprecision(10);
    ofs_v << setprecision(10);
    
    Wannierfile p2_WF;
    Wannierfile p_WF;
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
            if(cel.ReadWC(ifs_wan)==0)break;
            cel.organize_pos();
            cel.organize_wan();
            //read and write water_dipole
            Wannierfile WF;
            WF.Read_Wannier(cel);
            WF.Read_Dipole(cel,read_type);
            WF.Print_Dipole(ofs);
            //recort dt
            if(INPUT.ss_n++ == 1)
                INPUT.dt = WF.time - p_WF.time;
            //read and write water_vdipole and record tot_vdipole of each snapshot
            WF.Read_Vdipole(p_WF,p2_WF,ofs_v);
        }
    }
    
    cout << " dt: " << INPUT.dt << " total snapshot: "<< INPUT.ss_n << endl;
    INPUT.ss_n -= 2;
}

void Wannierfile::Routine_Ir_1(ofstream &ofs, ofstream &ofs_v, Vector3<double> *vdipole_file, int read_type)
{
    INPUT.ifs_cel.seekg(INPUT.cel_top);
    ifstream ifs(INPUT.geo_file.c_str());
    ifstream ifs_wan(INPUT.wan_file.c_str());
    file_assert(ifs,INPUT.geo_file);
    file_assert(ifs_wan,INPUT.wan_file);
    ofs << setprecision(10);
    ofs_v << setprecision(10);
    
    Wannierfile p2_WF;
    Wannierfile p_WF;
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
            if(cel.ReadWC(ifs_wan)==0)break;
            cel.organize_pos();
            cel.organize_wan();
            //read and write water_dipole
            Wannierfile WF;
            WF.Read_Wannier(cel);
            WF.Read_Dipole(cel,read_type);
            WF.Print_Dipole(ofs);
            //recort dt
            if(INPUT.ss_n++ == 1)
                INPUT.dt = WF.time - p_WF.time;
            //read and write water_vdipole and record tot_vdipole of each snapshot
            WF.Read_Vdipole(p_WF,p2_WF,ofs_v,vdipole_file[INPUT.ss_n-3]);
        }
    }
    cout << " dt: " << INPUT.dt << " total snapshot: "<< INPUT.ss_n << endl;
    INPUT.ss_n -= 2;
}

void TCF(Vector3<double> *A,double dt,int ss_n,int delta,int tjump,double t_cutoff,double alpha_cutoff)
{
    ofstream ofs("TCF.txt");
    ofs << setprecision(10);
    
    cout << "Calculating Time Correlation Function with time_cutoff = " << t_cutoff << ", alpha = " << alpha_cutoff << endl;
    for(int i=0; i<delta and i*(1+tjump)<ss_n-1;i++)
    {
        double t = dt*i*(1+tjump);
        if(t>t_cutoff)
        {
            double x=(t-t_cutoff)/(delta*dt*(1+tjump)-t_cutoff);
            //double fric=2*x*x*x-3*x*x+1;
            double fric = exp(-0.5*alpha_cutoff*alpha_cutoff*x*x);
            //cout << setw(10) << x << setw(10) <<fric << endl;
            ofs << t << " " << Time_correlation(A,ss_n-1,i*(1+tjump))*fric << endl;
        }
        else
            ofs << t << " " << Time_correlation(A,ss_n-1,i*(1+tjump)) << endl;
    }
}

int TCF2list(double *result, Vector3<double> *A,double dt,int ss_n,int delta,int tjump,double t_cutoff,double alpha_cutoff, bool output)
{
    if(output)
        cout << "Calculating Time Correlation Function with time_cutoff = " << t_cutoff << ", alpha = " << alpha_cutoff << endl;
    int i;
    for(i=0; i<delta and i*(1+tjump)<ss_n-1;i++)
    {
        double t = dt*i*(1+tjump);
        if(t>t_cutoff)
        {
            double x=(t-t_cutoff)/(delta*dt*(1+tjump)-t_cutoff);
            double fric = exp(-0.5*alpha_cutoff*alpha_cutoff*x*x);
            result[i]+=Time_correlation(A,ss_n-1,i*(1+tjump))*fric;
        }
        else
            result[i]+=Time_correlation(A,ss_n-1,i*(1+tjump));
    }
    return i;
}

int TCF2list2(double *result, Vector3<double> *A, Vector3<double> *B,double dt,int ss_n,int delta,int tjump,double t_cutoff,double alpha_cutoff,bool output )
{
    if(output)
        cout << "Calculating Time Correlation Function with time_cutoff = " << t_cutoff << ", alpha = " << alpha_cutoff << endl;
    int i;
    for(i=0; i<delta and i*(1+tjump)<ss_n-1;i++)
    {
        double t = dt*i*(1+tjump);
        if(t>t_cutoff)
        {
            double x=(t-t_cutoff)/(delta*dt*(1+tjump)-t_cutoff);
            double fric = exp(-0.5*alpha_cutoff*alpha_cutoff*x*x);
            result[i]+=Time_correlation(A,B,ss_n-1,i*(1+tjump))*fric;
        }
        else
            result[i]+=Time_correlation(A,B,ss_n-1,i*(1+tjump));
    }
    return i;
}

int TCF2list3(double *result, Vector3<double> *A, Vector3<double> *B, Vector3<double> *posi, Vector3<double> *posj, Vector3<double> *celldm, double pos_start, double pos_end, double dt,int ss_n,int delta,int tjump,double t_cutoff,double alpha_cutoff,bool output )
{
    assert(pos_start < pos_end);
    if(output)
        cout << "Calculating Time Correlation Function with time_cutoff = " << t_cutoff << ", alpha = " << alpha_cutoff << endl;
    int i;
    bool *condition; // to mark if the two fit our condition
    condition = new bool [ss_n];
    for(i=0; i<delta and i*(1+tjump)<ss_n-1;i++)
    {
        double t = dt*i*(1+tjump);
        //calculate conditions
        for(int j=0; j<ss_n - i*(1+tjump); j++)
        {
            double tmpdis = posi[j].distance_BC(posj[j+i*(1+tjump)],celldm[j]);
            if( tmpdis >= pos_start and tmpdis < pos_end )
                condition[j] = true;
            else
                condition[j] = false;
        }
        if(t>t_cutoff)
        {
            double x=(t-t_cutoff)/(delta*dt*(1+tjump)-t_cutoff);
            double fric = exp(-0.5*alpha_cutoff*alpha_cutoff*x*x);
            result[i]+=Time_correlation(A,B,condition,ss_n-1,i*(1+tjump))*fric;
        }
        else
            result[i]+=Time_correlation(A,B,condition,ss_n-1,i*(1+tjump));
    }
    delete [] condition;
    return i;
}

int TCF2list4(double **result, Vector3<double> *A, Vector3<double> *B, Vector3<double> *posi, Vector3<double> *posj, Vector3<double> *celldm, double pos_start, double pos_end, int pos_n, double dt,int ss_n,int delta,int tjump,double t_cutoff,double alpha_cutoff,bool output )
{
    assert(pos_start < pos_end);
    assert(pos_n > 0);
    if(output)
        cout << "Calculating Time Correlation Function with time_cutoff = " << t_cutoff << ", alpha = " << alpha_cutoff << endl;
    int i;
    double pos_step = (pos_end - pos_start)/pos_n;
    int *condition;
    bool *condition2; // to mark if the two fit our condition
    condition = new int [ss_n];
    condition2 = new bool [ss_n];
    for(i=0; i<delta and i*(1+tjump)<ss_n-1;i++)
    {
        double t = dt*i*(1+tjump);
        //calculate conditions
        for(int j=0; j<ss_n - i*(1+tjump); j++)
        {
            double tmpdis = posi[j].distance_BC(posj[j+i*(1+tjump)],celldm[j]);
            if( tmpdis >= pos_start and tmpdis < pos_end )
                condition[j] = (tmpdis - pos_start)/pos_step;
            else
                condition[j] = -1;
        }
        for(int j=0; j<pos_n; j++)
        {
            for(int k=0; k<ss_n - i*(1+tjump); k++)
            {
                if(condition[k]==j) condition2[k] = true;
                else condition2[k] = false;
            }
            if(t>t_cutoff)
            {
                double x=(t-t_cutoff)/(delta*dt*(1+tjump)-t_cutoff);
                double fric = exp(-0.5*alpha_cutoff*alpha_cutoff*x*x);
                result[j][i]+=Time_correlation(A,B,condition2,ss_n-1,i*(1+tjump))*fric;
            }
            else
                result[j][i]+=Time_correlation(A,B,condition2,ss_n-1,i*(1+tjump));
        }
    }
    delete [] condition;
    delete [] condition2;
    return i;
}

int TCF2list5(double *result, Vector3<double> *A, Vector3<double> *B, Vector3<double> *C, Vector3<double> *D,Vector3<double> *posi, Vector3<double> *posj, Vector3<double> *celldm, double pos_start, double pos_end, double dt,int ss_n,int delta,int tjump,double t_cutoff,double alpha_cutoff,bool output )
{
    assert(pos_start < pos_end);
    if(output)
        cout << "Calculating Time Correlation Function with time_cutoff = " << t_cutoff << ", alpha = " << alpha_cutoff << endl;
    int i;
    bool *condition; // to mark if the two fit our condition
    condition = new bool [ss_n];
    for(i=0; i<delta and i*(1+tjump)<ss_n-1;i++)
    {
        double t = dt*i*(1+tjump);
        //calculate conditions
        for(int j=0; j<ss_n - i*(1+tjump); j++)
        {
            double tmpdis = posi[j].distance_BC(posj[j+i*(1+tjump)],celldm[j]);
            if( tmpdis >= pos_start and tmpdis < pos_end and C[i]*D[j] >0)
                condition[j] = true;
            else
                condition[j] = false;
        }
        if(t>t_cutoff)
        {
            double x=(t-t_cutoff)/(delta*dt*(1+tjump)-t_cutoff);
            double fric = exp(-0.5*alpha_cutoff*alpha_cutoff*x*x);
            result[i]+=Time_correlation(A,B,condition,ss_n-1,i*(1+tjump))*fric;
        }
        else
            result[i]+=Time_correlation(A,B,condition,ss_n-1,i*(1+tjump));
    }
    delete [] condition;
    return i;
}

void FT(double *y, double dt, int n2, string out_file)
{
    ofstream ofs(out_file.c_str());
    assert(INPUT.T);
    if(INPUT.ensemble == "npt") Reorganize::Read_average_cel(1);
    assert(INPUT.Celldm.norm());
    
    //unit convert
    double Kb = 1.38064852*1e-23;  // boltzman constant J/K
    double C  = 299792458;      // speed of light m/s
    double Miu = 4*PI*1e-7;    // magnetic constant H/m or kg/s^2/A^2
    double D2C = 0.20819434*1.6*1e-19*1e-10;  // unit conv of dipole from debye to C*M
    double A2S = 2.418884326505*1e-17; // unit conv of time from a.u. time to s
    double P2S = 1e-12;     // unit conv of time from ps to s
    double M2C = 1e-5;      // unit conv of final output from m^-1 to 10^3 cm^-1
    
    
    //values in SI
    double SIT = INPUT.T;     // Temperature K
    double SIV = INPUT.Celldm.x*INPUT.Celldm.y*INPUT.Celldm.z*pow(INPUT.unitconv,3)* 1e-30; // Volume m^3
    
    //double uc = 36.0*PI/1.38/9.0 *0.20819434*0.20819434*1.6*1.6/2.418884326505/2.418884326505*pow(10.0,13)/INPUT.T/INPUT.Celldm.x/INPUT.Celldm.y/INPUT.Celldm.z/pow(INPUT.unitconv,3);
    double uc = Miu*C/(3*Kb*SIT*SIV)*D2C*D2C/A2S/A2S*P2S*M2C;
    
    double dk = 2*PI/dt/n2;
    double n3 = 3*dt*n2*INPUT.upper_limit/100;
    int thicken = INPUT.thicken;
    // calculate more points if we can by default (assuming the tcf(large t)=0, so that we can add more 0 in TCF)
    if(thicken==0)
    {
        if(n3 < n2)
            thicken = n2/n3;
        else
            thicken = 1;
    }
    dk = dk/thicken;
    
    cout << "Calculating Infrared Spectra with thicken = "  << thicken << endl;
    for(int i=0; i<=n2*thicken and i<=n3*thicken ; i++)
    {
        //dk*i/6/PI*100 unit convert from 2PI*ps-1 to cm-1
        ofs << dk*i/6/PI*100 << " " << uc*Fourier_transform(dk*i,y,dt,n2) << endl;
    }
}
