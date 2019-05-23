#include "Physics_structure_factor.h"


Structure_factor::Structure_factor(const Distributionfunction& in_D, double rou, double r0, int Ndr){
    auto G_r = in_D.get_y();
    auto _r = in_D.get_x();
    S_q.clear();
    _q.clear();
    compute_sturcture_factor(G_r.cbegin(),G_r.cend(),_r[1]-_r[0],std::back_inserter(S_q),rou,r0,Ndr);
    compute_q(std::back_inserter(_q));
}

void Structure_factor::Solve(g_citer begin, g_citer end, double dr, double rou, double r0, int Ndr){
    S_q.clear();
    _q.clear();
    compute_sturcture_factor(begin,end,dr,std::back_inserter(S_q),rou,r0,Ndr);
    compute_q(std::back_inserter(_q));
}


void Structure_factor::compute_sturcture_factor(g_citer begin, g_citer end, double dr, std::back_insert_iterator<std::vector<double>> result, double rou, double r0, int Ndr){
    
    // the number of point of both r and q
    if(Ndr == -1)
        N = end - begin;
    else
        N = Ndr;
    // dq in output
    dq = 2*PI/dr/N;
    
    double q=0;
    for(int j=0; j< N; ++j){
        double r=r0;
        double Sn=0;
        auto it=begin;
        for(int i=0; i< N and it!=end; ++i){
            Sn += r*sin(q*r)*(*it-1);
            r += dr;
            ++it;
        }
        *(result++)= (Sn==0) ? 1 : 4*PI*rou*Sn/q*dr+1;
        q += dq;
    }
}

void Structure_factor::compute_q(std::back_insert_iterator<std::vector<double>> result){
    assert(N>0);
    double q=0;
    for(int j=0; j< N; ++j){
        *(result++)=q;
        q+=dq;
    }
}
