#include "H_bondfile.h"
#include "input.h"

Water::Water()
{
}

Water::~Water()
{
}

Waterana::Waterana()
{
    allocate_pos = false;
    allocate_vnorm = false;
}

Waterana::Waterana(const Water& water)
{
    indexO = water.indexO;
    numH = water.numH;
    for (int i=0; i<numH; i++)
    {
        indexH[i] = water.indexH[i];
        disH[i] = water.disH[i];
    }
    allocate_pos = false;
    allocate_vnorm = false;
}

Waterana::Waterana(const Water& water,const Cellfile &cel,const double unitconv)
{
    indexO = water.indexO;
    numH = water.numH;
    for (int i=0; i<numH; i++)
    {
        indexH[i] = water.indexH[i];
        disH[i] = water.disH[i];
    }
    
    int O,H;
    for(int i=0; i<cel.ntype;i++) //find out index of O and H in Cel
    {
        if(cel.atoms[i].id == "O") O=i;
        if(cel.atoms[i].id == "H") H=i;
    }
    
    posO = cel.atoms[O].pos[indexO]*unitconv;
    for(int i=0; i<numH; i++)
    {
        posH.push_back(cel.atoms[H].pos[indexH[i]]*unitconv);
        OH.push_back(posH[i]-posO);
    }
    allocate_pos = true;
    set_OH();
    allocate_vnorm = false;
    
}

Waterana::Waterana(const Vector3<double> &inposO,const vector<Vector3<double> > &inposH)
{
    posO = inposO;
    posH = inposH;
    numH = posH.size();
    allocate_pos = true;
    set_OH();
    allocate_vnorm = false;
}



Waterana::Waterana(const Waterana& water)
{
    indexO = water.indexO;
    numH = water.numH;
    for (int i=0; i<numH; i++)
    {
        indexH[i] = water.indexH[i];
        disH[i] = water.disH[i];
    }
    
    allocate_pos = water.allocate_pos;
    allocate_vnorm = water.allocate_vnorm;
    posO = water.posO;
    posH = water.posH;
    OH = water.OH;
    vnorm = water.vnorm;
}

Waterana::~Waterana()
{
    
}

int Waterana::set_OH()
{
    assert(allocate_pos);
    if(OH.size() != numH)
    {
        OH.clear();
        for(int i=0; i<numH; i++)
            OH.push_back(posH[i]-posO);
        return numH;
    }
    else
        return 0;
    
}


int Waterana::set_vnorm(bool force_set)
{
    assert(numH==2);
    set_OH();
    if(!allocate_vnorm or force_set)
    {
        vnorm = OH[0].cross(OH[1]);
        vnorm /= vnorm.norm();
        allocate_vnorm = true;
        return 1;
    }
    else
        return 0;
}

Waterana Waterana::translate(Vector3<double> d, double mO, double mH)
{
    assert(allocate_pos);
    Vector3<double> tmpposO;
    vector<Vector3<double> > tmpposH;
    
    tmpposO = posO + d;
    for(Vector3<double> H : posH)
        tmpposH.push_back(H + d);
    
    Waterana tmp(tmpposO, tmpposH);
    return tmp;
}

Waterana Waterana::liberate(double d, double mO, double mH)
{
    assert(allocate_pos);
    if(numH != 2)
    {
        cerr << "Can not liberate a water with " << numH <<  " H atoms!" << endl;
        return *this;
    }
    
    set_vnorm();
    
    Vector3<double> vH1,vH2,vO,tmpposO;
    vector<Vector3<double> > tmpposH;
    
    vH1 = OH[0].cross(vnorm).uniform() * d;
    vH2 = OH[1].cross(vnorm).uniform() * d;
    vO = (vH1 * mH + vH2 * mH) / (mO * -1);
    
    tmpposO = posO + vO;
    tmpposH.push_back(posH[0] + vH1);
    tmpposH.push_back(posH[1] + vH2);
    Waterana tmp(tmpposO, tmpposH);
    return tmp;
}

Waterana Waterana::liberate(Vector3<double> axis, double d, double mO, double mH)
{
    assert(allocate_pos);
    
    Vector3<double> v,tmpposO;
    vector<Vector3<double> > tmpposH;
    
    tmpposO = posO;
    for(Vector3<double> H : posH)
    {
        v = (H-posO).cross(axis).uniform() * d;
        tmpposH.push_back(H + v);
        tmpposO -= v*mH/mO;
    }
    
    Waterana tmp(tmpposO, tmpposH);
    return tmp;
}

Waterana Waterana::bend(double d,double mO,double mH)
{
    assert(allocate_pos);
    if(numH != 2)
    {
        cerr << "Can not bend a water with " << numH <<  " H atoms!" << endl;
        return *this;
    }
    
    set_vnorm();
    
    Vector3<double> vH1,vH2,vO,tmpposO;
    vector<Vector3<double> > tmpposH;
    
    vH1 = OH[0].cross(vnorm).uniform() * d;
    vH2 = OH[1].cross(vnorm).uniform() * -1 * d;
    vO = (vH1 * mH + vH2 * mH) / (mO * -1);
    
    tmpposO = posO + vO;
    tmpposH.push_back(posH[0] + vH1);
    tmpposH.push_back(posH[1] + vH2);
    Waterana tmp(tmpposO, tmpposH);
    return tmp;
}

Waterana Waterana::nbend(int steps, double d, double mO, double mH, bool print, ostream &os)
{
    assert(steps >= 0);
    string label;
    
    convstring(steps,label);
    if(print)
        print_xyz(os,label + "d");
    
    if(steps==0)
    {
        return *this;
    }
    else
    {
        *this = bend(d,mO,mH);
        return nbend(steps-1,d,mO,mH,print,os);
    }

}

Waterana Waterana::astretch(double d,double mO,double mH)
{
    assert(allocate_pos);
    if(numH != 2)
    {
        cerr << "Can not stretch a water with " << numH <<  " H atoms!" << endl;
        return *this;
    }
    
    Vector3<double> vH1,vH2,vO,tmpposO;
    vector<Vector3<double> > tmpposH;
    
    vH1 = OH[0].uniform() * d;
    vH2 = OH[1].uniform() * -1 * d;
    vO = (vH1 * mH + vH2 * mH) / (mO * -1);
    
    tmpposO = posO + vO;
    tmpposH.push_back(posH[0] + vH1);
    tmpposH.push_back(posH[1] + vH2);
    Waterana tmp(tmpposO, tmpposH);
    return tmp;
}

Waterana Waterana::sstretch(double d,double mO,double mH)
{
    assert(allocate_pos);
    if(numH != 2)
    {
        cerr << "Can not stretch a water with " << numH <<  " H atoms!" << endl;
        return *this;
    }
    
    Vector3<double> vH1,vH2,vO,tmpposO;
    vector<Vector3<double> > tmpposH;
    
    vH1 = OH[0].uniform() * d;
    vH2 = OH[1].uniform() * d;
    vO = (vH1 * mH + vH2 * mH) / (mO * -1);
    
    tmpposO = posO + vO;
    tmpposH.push_back(posH[0] + vH1);
    tmpposH.push_back(posH[1] + vH2);
    Waterana tmp(tmpposO, tmpposH);
    return tmp;
}

Waterana& Waterana::operator= (const Waterana &water)
{
    indexO = water.indexO;
    numH = water.numH;
    for (int i=0; i<numH; i++)
    {
        indexH[i] = water.indexH[i];
        disH[i] = water.disH[i];
    }
    
    allocate_pos = water.allocate_pos;
    allocate_vnorm = water.allocate_vnorm;
    posO = water.posO;
    posH = water.posH;
    OH = water.OH;
    vnorm = water.vnorm;
    return *this;
}

void Waterana::set(const Waterana &water)
{
    indexO = water.indexO;
    numH = water.numH;
    for (int i=0; i<numH; i++)
    {
        indexH[i] = water.indexH[i];
        disH[i] = water.disH[i];
    }
    
    allocate_pos = water.allocate_pos;
    allocate_vnorm = water.allocate_vnorm;
    posO = water.posO;
    posH = water.posH;
    OH = water.OH;
    vnorm = water.vnorm;
}

void Waterana::print_xyz(ostream &os,string mark)
{
    assert(allocate_pos);
    os << numH+1 << endl << mark << endl;
    os << "O " << posO << endl;
    for(int i=0; i< numH; i++)
    {
        os << "H " << posH[i] << endl;
    }
}
