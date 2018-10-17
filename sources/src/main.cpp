#include "gheader.h"
#include "glog.h"
#include "gmpi.h"
#include "mconst.h"
#include "mfunc.h"
#include "punit.h"
#include "pconst.h"
#include "dcell.h"


int main(int argc,char** argv)
{
    GMPI->init(argc,argv);
    GLOG->init();

    //Inputfile in("INPUT");
    //cout << setprecision(10);
    //in >> INPUT;
    
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            //if(strncmp(argv[i],"count",5) == 0 or strncmp(argv[i],"-c",2) == 0)//{INPUT.calculation="count";}
            //else if(strncmp(argv[i],"help",4) == 0 or strncmp(argv[i],"-h",2) == 0)
            //{Help::Routine(INPUT.calculation); return 1;}
            //else if(strncmp(argv[i],"--cal=",6) == 0)
            //{string tmp = argv[i];INPUT.calculation=tmp.substr(6);}
            //else if(strncmp(argv[i],"--type=",7) == 0)
            //{string tmp = argv[i];convstring(tmp.substr(7),INPUT.type);}
            //else if(strncmp(argv[i],"-t",2) == 0)
            //{i++; string tmp = argv[i];convstring(tmp,INPUT.type);}
            //else if(strncmp(argv[i],"--delta=",8) == 0)
            //{string tmp = argv[i];convstring(tmp.substr(8),INPUT.delta);}
            //else if(strncmp(argv[i],"-d",2) == 0)
            //{string tmp = argv[++i];convstring(tmp,INPUT.type);}
        }
    }
    
    {
        Dcell test;
        //test.CLP << 1,0,0,0,2,0,0,0,3;
        test.POS.resize(2,3);
        test.POS = Eigen::MatrixXd::Random(2,3);
        
        GLOG->out_setprecision(10);
        //GLOG->out_stream() << "NCOR=" << GMPI->ncore() << endl;
        //GLOG->out_stream() << Legendre2(l_angs) << endl;
        GLOG->out_stream() << test.CLP << endl;
        GLOG->out_stream() << test.POS << endl;
        GLOG->out_stream() << test.POS.row(1)*test.CLP.inverse() << endl;
        //test.POS.row(1)
        //GLOG->stream() << "RANK=" << GMPI->rank() << endl;
        //if(INPUT.calculation=="test")
        //{}
        //else if(INPUT.calculation=="density") {Density::Routine();}
        //else if(INPUT.calculation=="print_water") {Waterfile::Print_water();}
        //else if(INPUT.calculation=="constrain_water") {Waterfile::Add_constrain();}
       // else if(INPUT.calculation=="animate_water") {Waterfile::Animate_water();}
        //else if(INPUT.calculation=="reorganize") {Reorganize::Routine();}
        //else if(INPUT.calculation=="hbs"){WahbondsT::Routine();}
        //else if(INPUT.calculation=="wannier" or INPUT.calculation=="ir"){Wannierfile::Routine_Ir();}
        //else if(INPUT.calculation=="wair"){WaIR::Routine();}
        //else if(INPUT.calculation=="dftcell"){DFTCell::Routine();}
        //else if(INPUT.calculation=="wannier"){Wannierfile::Routine();}
        //else if(INPUT.calculation=="rdf"){Rdf::Routine();}
        //else if(INPUT.calculation=="count"){Count::Routine();}
        //else if(INPUT.calculation=="hba"){Hba::Routine();}
        //else if(INPUT.calculation=="wfa"){Wfa::Routine();}
        //else if(INPUT.calculation=="wfat"){Wfat::Routine();}
    }

    GLOG->end();
    GMPI->end();
}
/*
 Outputfile* p;
 {
 Outputfile* of = new Outputfile("POSCAR");
 //Outputfile of("POSCAR");
 p = of;
 }
 cout << p->name_only();
 */
