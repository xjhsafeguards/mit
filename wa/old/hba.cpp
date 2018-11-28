#include "input.h"
#include "hba.h"

Hba::Hba(Hbondfile& HBF_in, Cellfile& Cel_in): HBF(HBF_in), Cel(Cel_in)
{
    //Hbondfile& HBF = HBF_in;
    //Cellfile& Cel = Cel_in;
    n_HB = 0;
}

Hba::~Hba()
{
    
}

Hbafile::Hbafile(int N_in): N(N_in)
{
    dt = INPUT.dt;
    n_jump = 1;
    rOH = new Vector3<double> [N];
    
    allocate_rOH_tcf = false;
}

Hbafile::~Hbafile()
{
    delete [] rOH;
    if(allocate_rOH_tcf)
        delete [] rOH_tcf;
}

void Hba::Routine()
{
    cout << "Analyzing Hbond using type: " << INPUT.type << endl;
    switch (INPUT.type){
        case 0:
            cout << "Calculation Rotational Correlation Time of the O-H" << endl;
            break;
        case 1:
            cout << "Calculation Rotational Correlation Time of the O-H" << endl;
            break;
        default:
            cout << "Wrong type" << endl;
            break;
    }
    
    
    if(INPUT.type == 0 or INPUT.type == 1)
    {
        ifstream ifs(INPUT.geo_file.c_str());
        file_assert(ifs,INPUT.geo_file.c_str());
        ofstream ofs1("rOH_t.txt");
        file_assert(ofs1,"rOH_t.txt");
        ofstream ofs2("rOH_tcf.txt");
        file_assert(ofs2,"rOH_tcf.txt");
        Hbafile HBAF((INPUT.ss_stop-INPUT.ss_start)/INPUT.ss_step+1);
        
        INPUT.ss_n = 0;
        double pt = 0;
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
                Hbondfile HF;
                HF.Read_Hbond(cel);
                Hba AF(HF,cel);
                // calculate average rOH and output the result
                AF.Read_rOH();
                AF.Print_rOH(ofs1);
                HBAF.rOH[INPUT.ss_n]=AF.rOH;
                // record dt
                if(INPUT.ss_n == 0)
                    pt = cel.time;
                if(INPUT.ss_n++ == 1)
                    HBAF.dt = cel.time - pt;
            }
        }
        HBAF.N = INPUT.ss_n;//total snapshot counted
        HBAF.Calculate_rOH_tcf(INPUT.delta, ofs2);
        HBAF.Calculate_t2(ofs2);
    }
}

void Hba::Read_rOH()
{
    assert(HBF.allocate_Hbonds);
    
    rOH.set(0,0,0);
    n_HB = 0;
    for(int i=0; i<HBF.nwater; i++)
    {
        Hbond& tmp = HBF.Hbonds[i];
        n_HB += tmp.donate_count;
        for(int j=0; j<tmp.donate_count ; j++)
        {
            
            rOH = rOH + Cel.atoms[HBF.H].pos[tmp.index_donate[j]].shortest(Cel.atoms[HBF.O].pos[HBF.waters[i].indexO],Cel.celldm);
        }
    }
    rOH /= n_HB;
}

void Hba::Print_rOH(ofstream& ofs)
{
    assert(ofs.good());
    ofs << setw(10) << Cel.time;
    rOH.print(ofs);
}

void Hbafile::Calculate_rOH_tcf(int n)
{
    assert(!allocate_rOH_tcf);
    assert(n_jump);
    int Ntmp = N/n_jump;
    if(n<Ntmp)
        n_rOH_tcf = n;
    else
        n_rOH_tcf = Ntmp;
    //allocate space for rOH_tcf
    rOH_tcf = new double [n_rOH_tcf];
    allocate_rOH_tcf = true;
    //calculate value for Time correlation function
    for(int i=0; i<n_rOH_tcf; i++)
    {
        rOH_tcf[i] = Hbafile::Time_correlation_Legendre2(rOH,N,i*(n_jump));
        if(i != 0)
            rOH_tcf[i] /= rOH_tcf[0];
    }
}

void Hbafile::Calculate_rOH_tcf(int n, ofstream& ofs)
// with output ofs
{
    assert(ofs.good());
    assert(!allocate_rOH_tcf);
    assert(n_jump);
    int Ntmp = N/n_jump;
    if(n<Ntmp)
        n_rOH_tcf = n;
    else
        n_rOH_tcf = Ntmp;
    //allocate space for rOH_tcf
    rOH_tcf = new double [n_rOH_tcf];
    allocate_rOH_tcf = true;
    //calculate value for Time correlation function
    for(int i=0; i<n_rOH_tcf; i++)
    {
        rOH_tcf[i] = Hbafile::Time_correlation_Legendre2(rOH,N,i*(n_jump));
        if(i != 0)
            rOH_tcf[i] /= rOH_tcf[0];
        ofs << rOH_tcf[i] << endl;
    }
}

void Hbafile::Calculate_t2(ofstream& ofs)
{
    assert(ofs.good());
    assert(allocate_rOH_tcf);
    assert(n_jump);
    my_assert(dt != -1,"Using more than one snapthot");
    
    cout << "Calculate rotational correlation time t2 with" << endl;
    cout << " dt = " << dt*n_jump << " ; n_step = " << n_rOH_tcf << endl;
    // calculate t2
    double tmp = rOH_tcf[0];
    rOH_tcf[0] = 1;
    t2 = Integration(rOH_tcf,dt*n_jump,n_rOH_tcf);
    rOH_tcf[0] = tmp;
    
    cout << " t2 = " << t2 << endl;
    
    ofs << endl;
    ofs << " dt = " << dt*n_jump << " ; n_step = " << n_rOH_tcf << endl;
    ofs << " t2 = " << t2 << endl;
}

double Hbafile::Time_correlation_Legendre2(Vector3<double> *y, int n_max, int n_step)
{
    assert(n_max>n_step);
    assert(n_step>=0);
    
    double *y1;
    y1 = new double [n_max-n_step+1];
    for(int i=0; i<(n_max-n_step+1); i++)
        y1[i]= Legendre2(y[i]*y[i+n_step]);
    return Integration(y1,1.0/(n_max-n_step),n_max-n_step);
}
