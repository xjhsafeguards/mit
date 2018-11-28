#include "input.h"

Input INPUT;

Input::Input()
{  
    type = 0;
    ensemble = "nvt";
    
    ss_start = 1;
    ss_stop = 1000000;
    ss_step = 1;
    ss_n = -1;
    dt = -1.0;
    
    in_file = "none";
    cel_file = "none";
    geo_file = "none";
    wan_file = "none";
    out_file = "none";
    
    ntype = 0;
    nband = 0;
    T = 0;
    Celldm.set(0,0,0);
    
    cutoff = 1000000;
    alpha = 0;
    unitconv = 0.529177;
    OH_distance = 1.26;
    OO_distance = 3.5;
    HOO_angle = 30;
    OW_distance = 0.8;
    HW_distance = 0.8;
    
    upper_limit = 5000;
    thicken = 0;
    delta = 100;
    eps = 0.01;
    
    displacement = 0;
    atom_num_tot = 0;
    pos_start[0] = 0;
    pos_start[1] = 0;
    pos_start_n = 1;
    
    tcf_max = 0;
}

Input::~Input()
{
    if(ntype)
    {
        delete [] atom_num;
        delete [] atom_mass;
    }
}

void Input::Read_input(string inputfile)
{
    ifstream ifs(inputfile.c_str());
    
    if(ifs.good())
    cout<< "Reading file " << inputfile << endl;
    else
    {
        cout<< "Something went wrong when reading " << inputfile << endl;
        exit(0);
    }
    
    string v;
    while(ifs.good())
    {
        ifs >> v;
        
        if(v== "calculation") Read_value(ifs,calculation);
        else if(v== "type") Read_value(ifs,type);
        else if(v== "ensemble") Read_value(ifs,ensemble);
        else if(v== "in_file") Read_value(ifs,in_file);
        else if(v== "cel_file") Read_value(ifs,cel_file);
        else if(v== "out_file") Read_value(ifs,out_file);
        else if(v== "geo_file") Read_value(ifs,geo_file);
        else if(v== "wan_file") Read_value(ifs,wan_file);
        else if(v== "ss_start") Read_value(ifs,ss_start);
        else if(v== "ss_stop") Read_value(ifs,ss_stop);
        else if(v== "ss_step") Read_value(ifs,ss_step);
        else if(v== "ss_n") Read_value(ifs,ss_n);
        else if(v== "dt") Read_value(ifs,dt);
        else if(v== "cutoff") Read_value(ifs,cutoff);
        else if(v== "alpha") Read_value(ifs,alpha);
        else if(v== "thicken") Read_value(ifs,thicken);
        else if(v== "upper_limit") Read_value(ifs,upper_limit);
        else if(v== "unitconv") Read_value(ifs,unitconv);
        else if(v== "delta") Read_value(ifs,delta);
        else if(v== "displacement") Read_value(ifs,displacement);
        else if(v== "eps") Read_value(ifs,eps);
        else if(v== "atom_num_tot") Read_value(ifs,atom_num_tot);
        else if(v== "pos_start") Read_values(ifs,pos_start,2);
        else if(v== "pos_start_n") Read_value(ifs,pos_start_n);
        else if(v== "tcf_max") Read_value(ifs,tcf_max);
        else if(v== "index") Read_line(ifs,index,8);
        else if(v== "ntype")
        {
            Read_value(ifs,ntype);
            atom_num = new int [ntype];
            atom_mass = new double [ntype];
            atom_name = new string [ntype];
        }
        else if(v== "nband") Read_value(ifs,nband);
        else if(v== "atom_num")
        {
            assert(ntype!=0);
            Read_values(ifs,atom_num,ntype);
        }
        else if(v== "atom_name")
        {
            assert(ntype!=0);
            Read_values(ifs,atom_name,ntype);
        }
        else if(v== "atom_mass")
        {
            assert(ntype!=0);
            Read_values(ifs,atom_mass,ntype);
        }
        else if(v== "T") Read_value(ifs,T);
        else if(v== "celldm1")
        {
            
            Read_value(ifs,Celldm.x);
            if(Celldm.y == 0)Celldm.y = Celldm.x;
            if(Celldm.z == 0)Celldm.z = Celldm.x;
        }
        else if(v== "celldm2") Read_value(ifs,Celldm.y);
        else if(v== "celldm3") Read_value(ifs,Celldm.z);
        else if(v== "OH_distance") Read_value(ifs,OH_distance);
        else if(v== "OO_distance") Read_value(ifs,OO_distance);
        else if(v== "HOO_angle") Read_value(ifs,HOO_angle);
        else if(v== "OW_distance") Read_value(ifs,OW_distance);
        else if(v== "HW_distance") Read_value(ifs,HW_distance);
        else
        {
            char ch = v[0];
            if(ch != '!' and ch != '#')
                cout << "Discard input: " << v << endl;
            ifs.ignore(200,'\n');
        }
        v = "#";
        
    }
    
    if(atom_num_tot == 0)
        for(int i=0; i<ntype; i++)
        {
            atom_num_tot += atom_num[i];
        }
    
    if(ensemble == "npt")
    {
        ifs_cel.open(cel_file.c_str(),ios::in);
        cel_top = ifs_cel.tellg();
    }
    
    ifs.close();
}


void Input::Print2screen(void)
{
    cout << "in_file: " << in_file << endl;
    cout << "atom_num: " << atom_num[0] << " " << atom_num[1] << endl;
    cout << "Unitconv "<< unitconv << endl;
    cout << "celldm: ";
    Celldm.print2screen();
    
}

int NCOR = 1;
int RANK = 0;
ofstream ofs_log;



