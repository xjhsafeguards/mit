#include "H_bond.h"

Hbond::Hbond()
{
    accept_count = 0;
    donate_count = 0;
    HB_count = 0;
}

Hbond::~Hbond()
{
}

void Hbond::Print(ofstream &ofs)
{
    ofs << accept_count << '\t' << donate_count << '\t';
    for(int j=0; j<8; j++)
    {
        if(j<accept_count)
            ofs << index_accept[j] << '\t';
        else
            ofs << "-1\t";
    }
    for(int j=0; j<4; j++)
    {
        if(j<donate_count)
            ofs << index_donate[j] << '\t';
        else
            ofs << "-1\t";
    }
    ofs << endl;
}

