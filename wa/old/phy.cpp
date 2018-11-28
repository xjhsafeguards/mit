#include "phy.h"


/* ****** _TCF_ ****** */

void Tcf::set_up_calculation()
{
    my_throw( n_jump != 0, "Can't generate tcf with n_jump = 0");
    my_throw( N > 1, "Can't generate tcf with less than 2 snapshots");
    
    int Ntmp = (N-2)/n_jump;
    if(n_tcf > Ntmp)
        n_tcf = Ntmp; // Should have some reminder
    
    tcf.clear();
}

void Tcf::Calculate_tcf()
{
    set_up_calculation();
    for(int i=0; i<= n_tcf; i++)
        tcf.push_back(Time_correlation(data1,N-1,i*(n_jump)));
    allocate_tcf = true;
}

void Tcf::Calculate_tcf_legendre(int n)
{
    set_up_calculation();
    for(int i=0; i<= n_tcf; i++)
        tcf.push_back(Time_correlation_Legendre(data1,N-1,i*(n_jump),n));
    allocate_tcf = true;
}

ostream& Tcf::Print_tcf(ostream &os) const
{
    assert(allocate_tcf);
    int count = 0;
    for( double value : tcf)
    {
        os << setw(10) << count++ * n_jump * dt <<  setw(15) << value << endl;
    }
    return os;
}


