#include "main.h"

int main(int argc,char** argv)
{
    
#ifdef __MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NCOR);
    MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
#endif
    
    clock_t start_time,end_time;
    start_time = clock();
    
    INPUT.Read_input();
    
    string rank_n;
    convstring(RANK,rank_n);
    ofs_log.open("Core" + rank_n + ".log");
    
    
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            if(strncmp(argv[i],"count",5) == 0 or strncmp(argv[i],"-c",2) == 0){INPUT.calculation="count";}
            else if(strncmp(argv[i],"help",4) == 0 or strncmp(argv[i],"-h",2) == 0)
            {Help::Routine(); return 1;}
            else if(strncmp(argv[i],"--cal=",6) == 0)
            {string tmp = argv[i];INPUT.calculation=tmp.substr(6);}
            else if(strncmp(argv[i],"--type=",7) == 0)
            {string tmp = argv[i];convstring(tmp.substr(7),INPUT.type);}
            else if(strncmp(argv[i],"-t",2) == 0)
            {i++; string tmp = argv[i];convstring(tmp,INPUT.type);}
            else if(strncmp(argv[i],"--delta=",8) == 0)
            {string tmp = argv[i];convstring(tmp.substr(8),INPUT.delta);}
            else if(strncmp(argv[i],"-d",2) == 0)
            {string tmp = argv[++i];convstring(tmp,INPUT.type);}
        }
    }
    
    {
        if(INPUT.calculation=="test")
        {
            cout << setprecision(15);

            system("ls -l");
        }
        else if(INPUT.calculation=="density") {Density::Routine();}
        else if(INPUT.calculation=="print_water") {Waterfile::Print_water();}
        else if(INPUT.calculation=="constrain_water") {Waterfile::Add_constrain();}
        else if(INPUT.calculation=="animate_water") {Waterfile::Animate_water();}
        else if(INPUT.calculation=="reorganize") {Reorganize::Routine();}
        else if(INPUT.calculation=="hbs"){Hbondfile::Routine();}
        //else if(INPUT.calculation=="wannier" or INPUT.calculation=="ir"){Wannierfile::Routine_Ir();}
        else if(INPUT.calculation=="ir"){Wannierfile::Routine_Ir();}
        else if(INPUT.calculation=="wannier"){Wannierfile::Routine();}
        else if(INPUT.calculation=="rdf"){Rdf::Routine();}
        else if(INPUT.calculation=="count"){Count::Routine();}
        else if(INPUT.calculation=="hba"){Hba::Routine();}
        else if(INPUT.calculation=="wfa"){Wfa::Routine();}
        else if(INPUT.calculation=="wfat"){Wfat::Routine();}
    }

    if(RANK == 0)
        cout << "Job Done!" << endl;
    
    ofs_log << " --------------- " << endl;
    ofs_log << "      Finish     " << endl;
    ofs_log << " --------------- " << endl;
    end_time = clock();
    ofs_log << "Totle Time : " << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
    time_t time_end = std::time(NULL);
    ofs_log << "End Time: " << ctime(&time_end) << endl;
    
#ifdef __MPI
    MPI_Finalize();
#endif
    
    //time_t time_end = std::time(NULL);
    //cout << " Now: " << ctime(&time_end) << endl;
}

/*
 Vector3<double> tmp1(3,4,2);
 Vector3<double> tmp2(3,7,5);
 Vector3<double> tmp3(6,7,8);
 Matrix3<double> tmp(tmp1,tmp2,tmp3);
 Matrix3<double> tmpm(tmp);
 double d =2 ;
 cout << tmpm.row(1) << endl;
 cout << tmpm[2] << endl;
 cout << tmp/d << endl;
 try{
 throw runtime_error("test");
 }
 catch (runtime_error err)
 {
 cout << err.what() << endl;
 }
 
 cout << tmp.transposition() << endl;
 cout << tmp.inverse() << endl;
 */
