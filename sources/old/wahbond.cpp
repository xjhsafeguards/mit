#include "wahbond.h"
double Wahbonds::OO_distance = 3.5;//to find H_bond(SI)
double Wahbonds::HOO_angle = 30;//to find H_bond(degree)

void WahbondsT::Routine(){
    Help::help("hbs",INPUT.type);
    switch(INPUT.type){
        case 0:{
            Outputfile of1("average_Hbs.txt");
            Wahbonds::write_average_hbonds_header(of1.stream());
            Outputfile of2("each_Hbs.txt");
            Wahbonds::write_each_hbonds_header(of2.stream());
            Outputfile of3("ion_Hbs.txt");
            Wahbonds::write_ions_hbonds_header(of3.stream());
            Outputfile of4("Hbonds.txt");
            Cell cel(INPUT);
            WahbondsT HBT;
            while(cel.next()){
                Wahbonds HB(cel);
                HB.set_hbonds();
                HBT.add(HB);
                HB.write_average_hbonds(of1.stream());
                HB.write_each_hbonds(of2.stream());
                HB.write_ions_hbonds(of3.stream());
            }
            HBT.write_ion(of4.stream());
            HBT.write_water(of4.stream());
        }
        break;
        case 1:{
            Outputfile of1("average_Hbs.txt");
            Wahbonds::write_average_hbonds_header(of1.stream());
            Outputfile of2("each_Hbs.txt");
            Wahbonds::write_each_hbonds_header(of2.stream());
            Outputfile of4("Hbonds.txt");
            Cell cel(INPUT);
            WahbondsT HBT;
            while(cel.next()){
                Wahbonds HB(cel);
                HB.set_hbonds();
                HBT.add_water(HB);
                HB.write_average_hbonds(of1.stream());
                HB.write_each_hbonds(of2.stream());
            }
            HBT.write_water(of4.stream());
        }
        break;
        case 2:{
            Outputfile of1("average_Hbs.txt");
            Wahbonds::write_average_hbonds_header(of1.stream());
            Outputfile of2("each_Hbs.txt");
            Wahbonds::write_each_hbonds_header(of2.stream());
            double z1 = INPUT.parameter.at(0), z2 = INPUT.parameter.at(1);
            Outputfile of4("Hbonds_z"+ to_string(z1) + to_string(z2) +".txt");
            Cell cel(INPUT);
            WahbondsT HBT;
            while(cel.next()){
                Wahbonds HB(cel);
                HB.set_hbonds();
                HBT.add_water(HB,z1,z2);
                HB.write_average_hbonds(of1.stream());
                HB.write_each_hbonds(of2.stream());
            }
            HBT.write_water(of4.stream());
        }
            break;
        default:
            break;
    }
}



Wahbonds::Wahbonds(Cell& in_cel, bool ifsave): Waters(in_cel,ifsave){
}
void Wahbonds::write_average_hbonds_header(ostream& os){
    os << setw(8) << "snapshot" << setw(15) << "time" << setw(15) << "avg_awater" << setw(15) << "avg_dwater" << setw(15) << "avg_n_water" << setw(15) << "avg_aion" << setw(15) << "avg_dion" << setw(15) << "avg_n_ion" << endl;
}
void Wahbonds::write_average_hbonds(ostream& os) const{
    os << setw(8) << cel->snapshot() << setw(15) << cel->time() << setw(15) << avg_awater << setw(15) << avg_dwater << setw(15) << avg_count_water << setw(15) << avg_aion << setw(15) << avg_dion << setw(15) << avg_count_ion << endl;
}
void Wahbonds::write_each_hbonds_header(ostream& os){
    os << setw(5) << "O";
    Hbond::write_detail_header(os);
}
void Wahbonds::write_each_hbonds(ostream& os) const{
    os << cel->snapshot() << '\t' << cel->time() << endl;
    for(int i=0; i<cel->anum("O"); i++){
        os << setw(5) << i;
        hbonds[i]->write_detail(os);
    }
}
void Wahbonds::write_ions_hbonds_header(ostream& os){
    os << setw(8) << "ss" << setw(15) << "time" << setw(5) << "O" << setw(5) << "numH";
    Hbond::write_detail_header(os);
}
void Wahbonds::write_ions_hbonds(ostream& os) const{
    int ss = cel->snapshot();
    double t = cel->time();
    for(int j=0;j<nion();j++){
        int index = get_ion(j);
        os << setw(8) << ss << setw(15) << t << setw(5) << index << setw(5) << get_nH(index);
        hbonds[index]->write_detail(os);
    }
}
void Wahbonds::init_hbonds(){
    if(!allocate_waters())
        set_waters();
    hbonds.clear();
    for(int i=0; i<nwater(); i++)
        hbonds.push_back(make_shared<Hbond>());
}
void Wahbonds::set_hbonds(){
    init_hbonds();
    int nO = nwater();
    int nn = nion();
    for(int i=0; i<nO; i++)
        for(int j=0; j<nO; j++)
            donate(i,j);
    if(nn != 0)
    {
        avg_aion = avg_aion/nn;
        avg_dion = avg_dion/nn;
        avg_count_ion = avg_dion + avg_aion;
    }
    else
    {
        avg_aion = -1;
        avg_dion = -1;
        avg_count_ion = -1;
    }
    avg_awater = avg_awater/(nO-nn);
    avg_dwater = avg_dwater/(nO-nn);
    avg_count_water = avg_awater + avg_dwater;
}
void Wahbonds::donate(int O1,int O2){
    if(O1==O2) return;
    int iO1 = get_indexO(O1);
    int iO2 = get_indexO(O2);
    if(cel->distance("O",iO1,"O",iO2) < OO_distance){
        for(int i=0; i<get_nH(O1); i++){
            double angle = cel->angle("H",get_indexH(O1,i),"O",iO1,"O",iO2);
            if( angle < HOO_angle){
                hbonds[O1]->push_back_angle(angle);
                hbonds[O1]->push_back_donate(O2);
                hbonds[O2]->push_back_accept(O1);
                if(waters[O1]->nH() == 2){
                    avg_dwater++;
                }
                else
                    avg_dion++;
                if(waters[O2]->nH() == 2)
                    avg_awater++;
                else
                    avg_aion++;
            }
        }
    }
}
void Wahbonds::test(bool condition,string information) const{
    if(!condition)
        throw(runtime_error("Wahbonds::Fail_"+information));
}
void WahbondsT::write_ion(ostream& os){
    os << " --- THE H-BOND NETWORK OF IONS --- " << endl;
    os << setw(10) << " HB_number" << setw(20) << "accept_count" << setw(20) << "percentage(%)" << endl;
    for(int i=0;i<10;i++)
        os << setw(10) << i << setw(20) << count_iona[i] << setw(20) << static_cast<double>(count_iona[i])/countion*100 << endl;
    os << setw(10) << " HB_number" << setw(20) << "donate_count" << setw(20) << "percentage(%)" << endl;
    for(int i=0;i<10;i++)
        os << setw(10) << i << setw(20) << count_iond[i] << setw(20) << static_cast<double>(count_iond[i])/countion*100 << endl;
    os << " Counting for ions = " << countion << endl;
    set_totaion();
    set_totdion();
    os << " Average accept HBs for ions is " << static_cast<double>(totaion)/countion << endl;
    os << " Average donate HBs for ions is " << static_cast<double>(totdion)/countion << endl;
}
void WahbondsT::write_water(ostream& os){
    os << " --- THE H-BOND NETWORK OF WATER MOLECLUES --- " << endl;
    os << setw(10) << " HB_number" << setw(20) << "HB_count" << setw(20) <<"percentage(%)" << endl;
    for(int i=0;i<10;i++)
        os << setw(10) << i << setw(20) << count_water[i] << setw(20) << static_cast<double>(count_water[i])/count*100 << endl;
    set_tothwater();
    os << " Average accept HBs for water molecules is " << static_cast<double>(tothwater-totaion+totdion)/2/count << endl;
    os << " Average donate HBs for water molecules is " << static_cast<double>(tothwater+totaion-totdion)/2/count << endl;
}
void WahbondsT::add(const Wahbonds& WF){
    add_water(WF);
    for(int j=0;j<WF.nion();j++){
        int index = WF.get_ion(j);
        count_water[WF.nhbond(index)]--;
        count_iona[WF.naccept(index)]++;
        count_iond[WF.ndonate(index)]++;
        countion++;
        count--;
    }
}
void WahbondsT::add_water(const Wahbonds& WF){
    count += WF.nwater();
    for(int j=0;j<WF.nwater();j++)
        count_water[WF.nhbond(j)]++;
}
void WahbondsT::add_water(const Wahbonds& WF, double z1, double z2){
    for(int j=0;j<WF.nwater();j++){
        double z = WF.posO_incell(j)[2];
        if(z>z1 and z<z2){
            count_water[WF.nhbond(j)]++;
            ++count;
        }
    }
}
void WahbondsT::set_totaion(){
    totaion = 0;
    for( auto it=count_iona.cbegin(); it!=count_iona.cend(); it++)
        totaion += it->first*it->second;
}
void WahbondsT::set_totdion(){
    totdion = 0;
    for( auto it=count_iond.cbegin(); it!=count_iond.cend(); it++)
        totdion += it->first*it->second;
}
void WahbondsT::set_tothwater(){
    tothwater = 0;
    for( auto it=count_water.cbegin(); it!=count_water.cend(); it++)
        tothwater += it->first*it->second;
}
