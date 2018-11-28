#include "wfa.h"
#include "input.h"

void Wfa::Routine()
{
    cout << " Analyze  water file using type: " << INPUT.type << endl;
    switch (INPUT.type) {
        case 0:
        default:
            cout << "Wrong type" << endl;
            break;
    }
}
void Wfat::Routine()
{
    cout << " Analyze time dependent water file using type: " << INPUT.type << endl;
    
    if(RANK == 0)
        switch (INPUT.type) {
            case 0:
                cout << "Calculation average 2nd order Rotational Correlation Time of the O-H covalent bond! " << endl;
                break;
            case 1:
                cout << "Calculation 2nd order Rotational Correlation Time of the O-H covalent bond of O: " << INPUT.index[0] << " H: " << INPUT.index[1] << endl;
                break;
            case 2:
                cout << "Calculation average 1st order Rotational Correlation Time of the O-H covalent bond! " << endl;
                break;
            case 3:
                cout << "Calculation average 2nd order Rotational Correlation Time of the H2O out of plane vector! " << endl;
                break;
            default:
                cout << "Wrong type" << endl;
                break;
        }
    
    switch (INPUT.type) {
        case 0:
        {
            int numO = INPUT.Atom_num("O");
/*
#ifdef __MPI
            int numOLoc = numO/NCOR;
            int numOLeft = numO%NCOR;
            int indexOInit = 0;
            if( RANK < numOLeft )
            {
                numOLoc += 1;
                indexOInit += RANK * numOLoc;
            }
            else if(RANK >= numO)
            {
                numOLoc = 0;
            }
            else
            {
                indexOInit = RANK * numOLoc + numOLeft
            }
#else
            
#endif */
            ifstream ifs(INPUT.geo_file.c_str());
            file_assert(ifs,INPUT.geo_file.c_str());
            ofstream ofs1("rOH_t.txt");
            file_assert(ofs1,"rOH_t.txt");
            ofs1 << setprecision(10);
            ofs1 << setw(5) << "O" << setw(15) << "time" << " rOH" << endl;
            ofstream ofs2("rOH_tcf.txt");
            file_assert(ofs2,"rOH_tcf.txt");
            ofs2 << setprecision(10);
            
            vector<Wfat> v_WFAT(numO,INPUT);
            
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
                    Waterfile WF;
                    WF.Read_water(cel);
                    Wfa AF(WF,cel,INPUT);
                    for(int j=0; j<WF.Show_nwater(); j++)
                    {
                        AF.Read_rOH(j);
                        Print_rOH(ofs1 << setw(5) << j << setw(15) << AF.Show_time() << " ", AF) << endl;
                        v_WFAT[j].Read_rOH(AF);
                    }
                }
            }
            Wfat& result = v_WFAT[0];
            result.Calculate_rOH_tcf();
            int c = 0;
            for( auto p = ++v_WFAT.begin(); p!=v_WFAT.end(); ++p)
            {
#ifndef __MPI
                cout << "Calculating tcf for atom: " << ++c << endl;
#endif
                (*p).Calculate_rOH_tcf();
                result.Sum_rOH_tcf(*p);
            }
            result.Div_rOH_tcf(numO);
            result.Print_rOH_tcf(ofs2);
            result.Calculate_t2();

            
        }
            break;
        case 1:
        {
            ifstream ifs(INPUT.geo_file.c_str());
            file_assert(ifs,INPUT.geo_file.c_str());
            ofstream ofs1("rOH_t.txt");
            file_assert(ofs1,"rOH_t.txt");
            ofs1 << setprecision(10);
            ofstream ofs2("rOH_tcf.txt");
            file_assert(ofs2,"rOH_tcf.txt");
            ofs2 << setprecision(10);
            
            Wfat WFAT(INPUT);
            
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
                    Waterfile WF;
                    WF.Read_water(cel);
                    //Wfa AF(WF,cel,INPUT);
                    Wfa AF(WF,cel,INPUT);
                    // calculate average rOH and output the result
                    AF.Read_rOH(INPUT.index[0]);
                    AF.Print_rOH(ofs1) << endl;
                    WFAT.Read_rOH(AF);
                }
            }
            WFAT.Calculate_rOH_tcf();
            WFAT.Print_rOH_tcf(ofs2);
            WFAT.Calculate_t2();
        }
        case 2:
        {
            int numO = INPUT.Atom_num("O");
            
            ifstream ifs(INPUT.geo_file.c_str());
            file_assert(ifs,INPUT.geo_file.c_str());
            ofstream ofs1("rOH_t.txt");
            file_assert(ofs1,"rOH_t.txt");
            ofs1 << setprecision(10);
            ofs1 << setw(5) << "O" << setw(15) << "time" << " rOH" << endl;
            ofstream ofs2("rOH_tcf_1.txt");
            file_assert(ofs2,"rOH_tcf_1.txt");
            ofs2 << setprecision(10);
            
            vector<Wfat> v_WFAT(numO,INPUT);
            
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
                    Waterfile WF;
                    WF.Read_water(cel);
                    Wfa AF(WF,cel,INPUT);
                    for(int j=0; j<WF.Show_nwater(); j++)
                    {
                        AF.Read_rOH(j);
                        Print_rOH(ofs1 << setw(5) << j << setw(15) << AF.Show_time() << " ", AF) << endl;
                        v_WFAT[j].Read_rOH(AF);
                    }
                }
            }
            Wfat& result = v_WFAT[0];
            result.Calculate_rOH_tcf_1();
            int c = 0;
            for( auto p = ++v_WFAT.begin(); p!=v_WFAT.end(); ++p)
            {
#ifndef __MPI
                cout << "Calculating tcf for atom: " << ++c << endl;
#endif
                (*p).Calculate_rOH_tcf_1();
                result.Sum_rOH_tcf(*p);
            }
            result.Div_rOH_tcf(numO);
            result.Print_rOH_tcf(ofs2);
            result.Calculate_t2();
            
            
        }
            break;
        case 3:{
            int numO = INPUT.Atom_num("O");
            ifstream ifs(INPUT.geo_file.c_str());
            file_assert(ifs,INPUT.geo_file.c_str());
            ofstream ofs1("rOH_ofp_t.txt");
            file_assert(ofs1,"rOH_ofp_t.txt");
            ofs1 << setprecision(10);
            ofs1 << setw(5) << "O" << setw(15) << "time" << " rOH" << endl;
            ofstream ofs2("rOH_ofp_tcf.txt");
            file_assert(ofs2,"rOH_ofp_tcf.txt");
            ofs2 << setprecision(10);
            
            vector<Wfat> v_WFAT(numO,INPUT);
            
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
                    Waterfile WF;
                    WF.Read_water(cel);
                    Wfa AF(WF,cel,INPUT);
                    for(int j=0; j<WF.Show_nwater(); j++)
                    {
                        AF.Read_rOH_ofp(j);
                        Print_rOH(ofs1 << setw(5) << j << setw(15) << AF.Show_time() << " ", AF) << endl;
                        v_WFAT[j].Read_rOH(AF);
                    }
                }
            }
            Wfat& result = v_WFAT[0];
            result.Calculate_rOH_tcf();
            int c = 0;
            for( auto p = ++v_WFAT.begin(); p!=v_WFAT.end(); ++p)
            {
#ifndef __MPI
                cout << "Calculating tcf for atom: " << ++c << endl;
#endif
                (*p).Calculate_rOH_tcf();
                result.Sum_rOH_tcf(*p);
            }
            result.Div_rOH_tcf(numO);
            result.Print_rOH_tcf(ofs2);
            result.Calculate_t2();
        }
            break;
    }
}


Wfa::Wfa(const Waterfile& WF_in, const Cellfile& Cel_in): WF(WF_in), Cel(Cel_in)
{
    n_OH = 0;
    unitconv = 1;
}

Wfa::Wfa(const Waterfile& WF_in, const Cellfile& Cel_in, const Input& INPUT_in): WF(WF_in), Cel(Cel_in)
{
    n_OH = 0;
    unitconv = INPUT_in.unitconv;
}

Wfa::~Wfa()
{
    
}

void Wfa::Read_rOH()
{
    assert(WF.allocate_waters);
    
    rOH.set(0,0,0);
    n_OH = 0;
    for(int i=0; i<WF.nwater; i++)
    {
        n_OH += WF.waters[i].numH;
        for(int j=0; j<WF.waters[i].numH ; j++)
        {
            rOH += WF.posH(i,j).shortest(WF.posO(i),Cel.celldm)*unitconv;
        }
    }
    rOH /= n_OH;
}

void Wfa::Read_rOH(int index_O)
{
    assert(WF.allocate_waters);
    assert(index_O < WF.nwater);
    
    rOH.set(0,0,0);
    n_OH = WF.Show_nH(index_O);
    //cout << index_O <<" "<< n_OH << " " << unitconv << endl;
    for(int j=0; j<WF.waters[index_O].numH ; j++)
    {
        rOH += WF.posH(index_O,j).shortest(WF.posO(index_O),Cel.celldm)*unitconv;
    }
    rOH /= n_OH;
}

void Wfa::Read_rOH_ofp(int index_O)
{
    assert(WF.allocate_waters);
    assert(index_O < WF.nwater);
    //assert(WF.waters[index_O].numH == 2);
    if(WF.waters[index_O].numH == 2){
        rOH.set(0,0,0);
        n_OH = 1;
        rOH += (WF.posH(index_O,0).shortest(WF.posO(index_O),Cel.celldm)*unitconv).cross((WF.posH(index_O,1).shortest(WF.posO(index_O),Cel.celldm)*unitconv));
    }
    else if(WF.waters[index_O].numH > 2){
        rOH.set(0,0,0);
        n_OH = 1;
        rOH += (WF.posH(index_O,0).shortest(WF.posO(index_O),Cel.celldm)*unitconv).cross((WF.posH(index_O,1).shortest(WF.posO(index_O),Cel.celldm)*unitconv));
        cout << "Warning more than two H detected while Read_rOH_ofp" << endl;
    }
    else if(WF.waters[index_O].numH < 2){
        rOH.set(0,0,0);
        n_OH = 1;
        cout << "Warning less than two H detected while Read_rOH_ofp" << endl;
    }
    
}

void Wfa::Read_rOH(int index_O,int index_H)
{
    assert(WF.allocate_waters);
    assert(index_O < WF.nwater);
    assert(index_H < WF.waters[index_O].numH);
    
    rOH.set(0,0,0);
    n_OH = 1;
    
    rOH = WF.posH(index_O,index_H).shortest(WF.posO(index_O),Cel.celldm)*unitconv;
}

Wfat::Wfat(): N(0), dt(-1), lt(0), allocate_rOH_tcf(false), n_rOH_tcf(0), n_jump(1) {}

Wfat::Wfat(const Input& INPUT_in) : Wfat()
{
    n_rOH_tcf = INPUT_in.tcf_max;
}

Wfat::Wfat(const Wfat& WFAT) :
N(WFAT.N), dt(WFAT.dt), lt(WFAT.lt), allocate_rOH_tcf(WFAT.allocate_rOH_tcf), n_rOH_tcf(WFAT.n_rOH_tcf), n_jump(WFAT.n_jump),
rOH_tcf(WFAT.rOH_tcf), t2(WFAT.t2) {}

Wfat& Wfat::operator=(const Wfat& WFAT)
{
    N=WFAT.N;
    dt=WFAT.dt;
    lt=WFAT.lt;
    allocate_rOH_tcf=WFAT.allocate_rOH_tcf;
    n_rOH_tcf=WFAT.n_rOH_tcf;
    n_jump=WFAT.n_jump;
    rOH_tcf=WFAT.rOH_tcf;
    t2=WFAT.t2;
    return *this;
}

Wfat::~Wfat()
{
    
}

Wfat& Wfat::operator<<(const Wfa &wfa)
{
    
    if(++N > 3)
        my_assert( abs(wfa.Cel.time - lt - dt) < 0.01, "Time interval between snapshots should be equal");
    else if(N == 2)
        dt = wfa.Cel.time - lt;
    //changed 2018.4.2 by jianhang
    rOH.push_back(wfa.rOH.uniform());
    lt = wfa.Cel.time;
    return *this;
}

Wfat& Wfat::Read_rOH(const Wfa &wfa)
{
    
    if(++N > 3)
        my_assert( abs(wfa.Cel.time - lt - dt) < 0.01, "Time interval between snapshots should be equal");
    else if(N == 2)
        dt = wfa.Cel.time - lt;
    //changed 2018.4.2 by jianhang
    rOH.push_back(wfa.rOH.uniform());
    lt = wfa.Cel.time;
    return *this;
}

void Wfat::Calculate_rOH_tcf()
{
    assert(n_jump);
    assert(N > 1);
    
    int Ntmp = (N-2)/n_jump; // max steps can perform, i.e. i*(n_jump) < N - 1
    if(n_rOH_tcf > Ntmp)
        n_rOH_tcf = Ntmp;
    
    rOH_tcf.clear();
    
    for(int i=0; i<= n_rOH_tcf; i++)
    {
        rOH_tcf.push_back(Time_correlation_Legendre2(rOH,N-1,i*(n_jump)));
        //if( i!= 0)
           //rOH_tcf[i] /= rOH_tcf[0];
    }
    //unify lost information
    //rOH_tcf[0] = 1;
    allocate_rOH_tcf = true;
}

void Wfat::Calculate_rOH_tcf(int delta)
{
    assert(n_jump);
    assert(N > 1);
    
    int Ntmp = (N-2)/n_jump; // max steps can perform, i.e. i*(n_jump) < N - 1
    (delta<Ntmp) ? n_rOH_tcf= delta : n_rOH_tcf = Ntmp;
    
    rOH_tcf.clear();
    
    for(int i=0; i<= n_rOH_tcf; i++)
    {
        rOH_tcf.push_back(Time_correlation_Legendre2(rOH,N-1,i*(n_jump)));
        //if( i!= 0)
            //rOH_tcf[i] /= rOH_tcf[0];
    }
    //unify lost information
    //rOH_tcf[0] = 1;
    allocate_rOH_tcf = true;
}

void Wfat::Calculate_rOH_tcf_1()
{
    assert(n_jump);
    assert(N > 1);
    
    int Ntmp = (N-2)/n_jump; // max steps can perform, i.e. i*(n_jump) < N - 1
    if(n_rOH_tcf > Ntmp)
        n_rOH_tcf = Ntmp;
    
    rOH_tcf.clear();
    correlate_function(rOH.begin(),rOH.end(),back_inserter(rOH_tcf),n_rOH_tcf+1,n_jump);
    allocate_rOH_tcf = true;
}

void Wfat::Calculate_t2()
{
    assert(allocate_rOH_tcf);
    assert(n_jump);
    
    my_assert(n_rOH_tcf > 1, " Cant calculate t2 with one tcf points");
    
    cout << "Calculate rotational correlation time t2 with" << endl;
    cout << " dt = " << dt*n_jump << " ; n_step = " << n_rOH_tcf << endl;
    // calculate t2
    t2 = Integration(rOH_tcf,dt*n_jump,n_rOH_tcf);
    
    cout << " t2 = " << t2 << endl;
    //ofs << endl;
    //ofs << " dt = " << dt*n_jump << " ; n_step = " << n_rOH_tcf << endl;
    //ofs << " t2 = " << t2 << endl;
}

ostream& Wfat::Print_rOH_tcf(ostream &os) const
{
    assert(allocate_rOH_tcf);
    int count = 0;
    for( double tcf : rOH_tcf)
    {
        os << setw(15) << count++ * n_jump * dt <<  setw(20) << tcf << endl;
    }
    return os;
}

Wfat& Wfat::Sum_rOH_tcf(const Wfat& WFAT)
{
    assert(allocate_rOH_tcf == WFAT.allocate_rOH_tcf);
    assert(n_rOH_tcf == WFAT.n_rOH_tcf );
    assert(n_jump == WFAT.n_jump);
    auto d2 = WFAT.rOH_tcf.begin();
    for(auto d1 = rOH_tcf.begin(); d1 != rOH_tcf.end(); ++d1, ++d2)
        *d1 += *d2;
    return *this;
}

Wfat& Wfat::Div_rOH_tcf(double d)
{
    assert(allocate_rOH_tcf);

    for(auto d1 = rOH_tcf.begin(); d1 != rOH_tcf.end(); ++d1)
        *d1 /= d;
    return *this;
}

/* bp
 N = 0;
 dt = -1;
 lt = 0;
 
 allocate_rOH_tcf = false;
 n_rOH_tcf = 0;
 n_jump = 1;
 */

