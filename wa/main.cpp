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
    
    string rank_n;
    convstring(RANK,rank_n);
    ofs_log.open("Core" + rank_n + ".log");

    Inputfile in("INPUT");
    cout << setprecision(10);
    in >> INPUT;
    
    if(argc != 1)
    {
        for(int i=1; i<argc; ++i)
        {
            if(strncmp(argv[i],"count",5) == 0 or strncmp(argv[i],"-c",2) == 0){INPUT.calculation="count";}
            else if(strncmp(argv[i],"help",4) == 0 or strncmp(argv[i],"-h",2) == 0)
            {Help::Routine(INPUT.calculation); return 1;}
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
            Outputfile* p;
            {
                Outputfile* of = new Outputfile("POSCAR");
                //Outputfile of("POSCAR");
                p = of;
            }
            cout << p->name_only();

        }
        //else if(INPUT.calculation=="density") {Density::Routine();}
        //else if(INPUT.calculation=="print_water") {Waterfile::Print_water();}
        //else if(INPUT.calculation=="constrain_water") {Waterfile::Add_constrain();}
       // else if(INPUT.calculation=="animate_water") {Waterfile::Animate_water();}
        //else if(INPUT.calculation=="reorganize") {Reorganize::Routine();}
        else if(INPUT.calculation=="hbs"){WahbondsT::Routine();}
        //else if(INPUT.calculation=="wannier" or INPUT.calculation=="ir"){Wannierfile::Routine_Ir();}
        else if(INPUT.calculation=="wair"){WaIR::Routine();}
        else if(INPUT.calculation=="dftcell"){DFTCell::Routine();}
        //else if(INPUT.calculation=="wannier"){Wannierfile::Routine();}
        //else if(INPUT.calculation=="rdf"){Rdf::Routine();}
        //else if(INPUT.calculation=="count"){Count::Routine();}
        //else if(INPUT.calculation=="hba"){Hba::Routine();}
        //else if(INPUT.calculation=="wfa"){Wfa::Routine();}
        //else if(INPUT.calculation=="wfat"){Wfat::Routine();}
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
    
}


/* test for gmath
 int i1=2;
 double d1 = 2;
 Vector3<double > x(2,0,0),y(0,2,0),z(0,0,2),t;
 Matrix3<double > m1(1,0,0,0,1,0,0,0,1),m2(x,y,z),m3(m1),m4;
 cout << x << endl;
 cout << "norm" << x.norm() << " " << y.norm() << " " << z.norm() << endl;
 cout << "unitform" << x.uniform() << endl;
 cout << "angle" << x.vangle(y) << endl;
 m4.set(1,2,3,4,5,5,6,7,8);
 cout << m1 << endl;
 cout << m2 << endl;
 cout << m3 << endl;
 cout << m4 << endl;
 cout << "row " <<m4.row(1) << endl;
 cout << "m*x x*m " << m4*x << " " << x*m4 << endl;
 cout << "-\n" << m1-m2 << endl;
 
 vector<double> te;
 for(int i=0; i<100; i++)
 te.push_back(i);
 for(int i=0; i<100; i++)
 cout << i << " "<< (correlate(te.cbegin(),te.cend(),i) == correlate2(te.cbegin(),te.cend(),te.cbegin(),te.cend(),i)) << endl;
 for(int i=0; i<100; i++)
 cout << i << " "<< (correlate_f(te.cbegin(),te.cend(),[](double d,double d2){return d*d2;},i) == correlate2_f(te.cbegin(),te.cend(),te.cbegin(),te.cend(),[](double d,double d2){return d*d2;},i)) << endl;
 
 vector<double > r1,r2;
 correlate_function(te.cbegin(),te.cend(),back_inserter(r1),100);
 correlate2_function(te.cbegin(),te.cend(),te.cbegin(),te.cend(),back_inserter(r2),100);
 correlate_f_function(te.cbegin(),te.cend(),back_inserter(r1),[](double d,double d2){return d*d2;},100);
 correlate2_f_function(te.cbegin(),te.cend(),te.cbegin(),te.cend(),back_inserter(r2),[](double d,double d2){return d*d2;},100);
 for(int i=0; i<200; i++)
 cout << i << " "<< (r1.at(i) == r2.at(i)) << endl;
 */

/* test for gfile
 Inputfile in("README.md"),in2;
 
 in.open("main.cpp");
 string ss,ss2;
 getline(in.stream(),ss);
 
 while (in)
 {
 //getline(in.stream(),ss);
 in >> ss >> ss2;
 cout << ss << " " << ss2 <<  endl;
 }
 
 cout << in.type() << endl;
 */

/* test for catom ccellp cwannier
 Inputfile in("example/water.pos"), in2("example/water.wfc");
 string buffer;
 cout << setprecision(10);
 
 //Atom at2("set");
 //cout << "name: " << at2.get_name();
 
 //in = in2;
 Atom at("O",64),at2;
 Wannier wan(256);
 at2 = at;
 at.set_label(in.stream());
 wan.set_label(in2.stream());
 cout << "time " << at.label(1) << endl;
 cout << "time " << wan.label(1) << endl;
 cout << at.num() << setw(10) << at.get_na() << endl;
 in >> at;
 in2 >> wan;
 cout << at.num() << setw(10) << at.get_na() << endl;
 cout << at[0]/Unit::Bohr2A << endl;
 for(auto it=at2.cbegin(); it!= at2.cend(); it++)
 cout << *it/Unit::Bohr2A << endl;
 cout << wan[0]/Unit::Bohr2A << endl;
 cout << at.isfractional() << endl;
 cout << wan.isfractional() << endl;
 //at.check_label(in.stream());
 //Cellp c;
 //in >> c;
 //cout << "test " << buffer;
 //cout << c[0] << c[1] << c[2] << endl;
 */

/* test for input
Inputfile in("INPUT");
cout << setprecision(10);
in >> INPUT;
cout << INPUT.cellp() << endl;
cout << INPUT.get_atom_num("O")<< endl;
cout << INPUT.get_atom_index("H")<< endl;
*/

/* test for cell
 cout << setprecision(15);
 Cell cel(INPUT);
 cel.skip();
 cel.set();
 cout << cel.anum() << endl;
 cout << cel.anum("H") << endl;
 cout << cel.aindex("O") << endl;
 cout << cel.apos("O",1)/Unit::Bohr2A << endl;
 cel.set_atoms_fractional();
 cout << cel.apos("O",1) << endl;
 cel.set_atoms_fractional(false);
 cout << cel.apos("O",1)/Unit::Bohr2A << endl;
 cout << cel.wpos(0)/Unit::Bohr2A  << endl;
 
 cout << "snapshot" << endl;
 for (auto it = cel.acbegin("O"); it!=cel.acend("O"); it++)
 cout << *it/Unit::Bohr2A << endl;
 
 cout << "cellp" << cel.box() << endl;
 
 /////////////////////////
 Cell cel(INPUT);
 cel.init();
 double uc = Unit::Bohr2A;
 while(cel.inrange()){
 cout << "snapshot" << cel.get_num() << endl;
 //for (auto it = cel.acbegin("O"); it!=cel.acend("O"); it++)
 //cout << *it/uc << endl;
 int i=0,j=1,k=1,l=1;
 cout << "1: " << cel.apos(i,j) << endl;
 cout << "1: " << cel.apos("O",j) << endl;
 cout << "f1:" << cel.box()*cel.afpos(i,j) << endl;
 cout << "2: " << cel.apos(k,l) << endl;
 cout << "2: " << cel.apos("H",l) << endl;
 cout << "f2: " << cel.box()*cel.afpos(k,l) << endl;
 cout << "distance1: "<< (cel.apos("O",j) - cel.apos("H",l)).norm() << endl;
 cout << "distance2: "<< cel.distance(i,j,k,l) << endl;
 //cout << "cellp" << endl;
 //cout << cel.box()/uc << endl;
 cel.next();
 //system("ls -l");
 */

/* test for wawater anwater
 Outputfile of("example/test.txt");
 Cell cel(INPUT);
 double uc = Unit::Bohr2A;
 while(cel.next()){
 Waters wa(cel);
 wa.set_waters();
 of << wa;
 }
 */

/* test for anhbond
 Hbond hb;
 hb.push_back_accept(11);
 hb.write_detail(cout);
 cout << hb.check_accept(10) << endl;
 cout << hb.check_accept(11) << endl;
 cout << hb.check_donate(11) << endl;
 cout << hb.check_bond(11) << endl;
 */

/* test for andipole wadipole
 //Vector3<double> x(14,0,0),y(0,23,0);
 //cout << acos(x*y/x.norm()/y.norm())/PI*180;
 Cell cel(INPUT);
 WaIR WDT;
 while(cel.next()){
 Wadipoles WD(cel);
 //WD.set_dipoles_wannier();
 //WD.write_each_wannier(cout);
 WD.set_dipoles_dipole();
 //WD.write_each_dipole(cout);
 WDT.add(WD);
 WDT.set_V(cel.volume());
 }
 //WDT.calculate_vdipole();
 //WDT.write_vdipole(cout);
 WDT.calculate_tcf(INPUT.delta);
 WDT.calculate_ir(INPUT.upper_limit);
 //WDT.write_tcf(cout);
 WDT.write_ir(cout);
  ////more
 cout << setprecision(10);
 Cell cel(INPUT);
 WadipolesT WDT;
 double VL;
 while(cel.next()){
 Wadipoles WD(cel);
 //WD.set_dipoles_wannier();
 //WD.write_each_wannier(cout);
 WD.set_dipoles_dipole();
 //WD.write_each_dipole(cout);
 WDT.add(WD);
 VL = cel.volume();
 cel.write_in(cout);
 }
 WDT.calculate_vdipole();
 auto vd = WDT.get_vdipole();
 IR ir1;
 IR3<double> ir;
 ir.set_dt(WDT.get_dt());
 ir1.set_dt(WDT.get_dt());
 ir.set_V(VL);
 ir1.set_V(VL);
 for( auto& it: *vd)
 cout << it << endl;
 ir1.set_vdipole(vd);
 ir1.calculate_tcf();
 ir.push_back_vdipole(vd);
 ir.push_back_condition(make_shared<vector<double> >(10));
 ir.calculate_tcf([](double a1, double a2){return 1;});
 ir.write_tcf(cout);
 ir1.write_tcf(cout);
 */
/* test for c* read cif
 Inputfile inf("example/0001.cif");
 Cellp cp;
 inf >> cp;
 //cp.set_parameter(6.74748725,8.14451453,10.49178793,90.00214572,90.00083977,96.88413802);
 cp.write_POSCAR(cout);
 double a,b,c,d,e,f;
 cp.get_parameter(a,b,c,d,e,f);
 cout << a << '\n' << b << '\n' << c << '\n' << d << '\n' << e << '\n' << f << endl;
 Atom at;
 at.set_na(40);
 at.set_name("H");
 inf >> at;
 at.write_in(cout);
 Cell cel;
 //DFTCell cel;
 inf >> cel;
 //cel.write_in(cout);
 Outputfile of("POSCAR");
 of << cel;
*/
