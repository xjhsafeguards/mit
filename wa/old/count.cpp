#include "count.h"
#include "input.h"

void Count::Routine()
{
    cout << "Counting the System with method: " << INPUT.type << endl;
    switch (INPUT.type){
        case 0:
        default:
        cout << "Reading snapshots of the system" << endl;
    }
    Count::Count_ss(INPUT.cel_file,4);
    Count::Count_ss(INPUT.geo_file,INPUT.atom_num_tot + 1);
    Count::Count_ss(INPUT.wan_file,INPUT.nband + 1);
}

void Count::Count_ss(string file, int lps)
{
    if(file != "none")
    {
        ifstream ifs(file.c_str());
        file_assert(ifs,file);
        int count_ss=0;
        
        while(ifs.good())
        {
            ifs.ignore(500,'\n');
            count_ss++;
        }
        
        count_ss /= lps;
        
        cout << file << " contains " << count_ss << " snapshots" << endl;
    }
}
