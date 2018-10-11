#include "wadipole.h"
double Wadipoles::OW_distance = 0.8;

Wadipoles::Wadipoles(Cell& in_cel,bool ifsave): Wahbonds(in_cel,ifsave){
}
void Wadipoles::write_each_wannier(ostream& os) const{
    int n = nwater();
    test(allocate_dipoles(n-1), "write_each_wannier::no_dipoles_allocated");
    for(int i=0; i<n; i++){
        os << setw(5) << i;
        dipoles.at(i)->write_wannier(os);
        os << endl;
    }
}
void Wadipoles::write_each_dipole(ostream& os) const{
    int n = nwater();
    test(allocate_dipoles(n-1), "write_each_dipole::no_dipoles_allocated");
    for(int i=0; i<n; i++){
        os << setw(5) << i << "  ";
        dipoles.at(i)->write_dipole(os);
        os << endl;
    }
}
void Wadipoles::init_dipoles(){
    if(!allocate_waters())
        set_waters();
    dipoles.clear();
    for(int i=0; i<nwater(); i++)
        dipoles.push_back(make_shared<Dipole>());
}
void Wadipoles::set_dipoles_wannier(){
    if(!allocate_dipoles())
        init_dipoles();
    for(int i=0; i<nwater(); i++){
        dipoles.at(i)->set_wannier(cel,"O",get_indexO(i),OW_distance);
        test(dipoles.at(i)->n()==4, "set_dipoles_wannier::water_No. " + to_string(i));
    }
}
void Wadipoles::set_dipoles_dipole(int type){
    if(!allocate_dipoles())
        set_dipoles_wannier();
    switch (type) {
        case 0:
            set_dipoles_total_dipole();
            break;
        default:
            break;
    }
}
void Wadipoles::test(bool condition,string information) const{
    if(!condition)
        throw(runtime_error("Wadipoles::Fail_"+information));
}
void Wadipoles::set_dipoles_total_dipole(){
    tot_dipole.set(0,0,0);
    for(int i=0; i<nwater(); ++i){
        if(get_nH(i) == 2){
            position_type tmpd;
            int O = get_indexO(i);
            for(int j=0; j<get_nH(i); ++j){
                tmpd += cel->cvector("O",O,"H",get_indexH(i,j));
            }
            for(int j=0; j<nwannier(i); ++j){
                tmpd -= cel->cvector("O",O,get_wannier(i,j))*2;
            }
        tot_dipole += tmpd;
        set_dipole(i,tmpd,0);
        }
    }
}


// WadipolesT
void WadipolesT::write_vdipole(ostream& os) const{
    os << setw(15) << dt  << setw(10) << vdipole.size() << endl ;
    for(const auto& v : vdipole)
        os << v << endl;
}

void WadipolesT::add(const Wadipoles& WD){
    set_time(WD.time());
    dipole.push_back(WD.get_tot_dipole());
}
void WadipolesT::add(const position_type& D, double in_time){
    set_time(in_time);
    dipole.push_back(D);
}
void WadipolesT::read_vdipole(istream& is){
    int size;
    position_type tmp;
    vdipole.clear();
    is >> dt >> size;
    for(int i=0;i<size;i++){
        is >> tmp;
        vdipole.push_back(tmp);
    }
}
void WadipolesT::calculate_vdipole(){
    test(allocate_dt(),"calculate_vdipole::no_dt");
    test(allocate_dipole(),"calculate_vdipole::no_dipole");
    derivate(dipole.cbegin(),dipole.cend(),back_inserter(vdipole),dt);
}
void WadipolesT::set_time(double nt){
    test(time==-1 or dt==-1 or abs(nt-time-dt) < 0.0000001, "set_time::different_time_interval::time" + to_string(nt));
    if(time != -1)
        dt = nt - time;
    time = nt;
}

void WadipolesT::test(bool condition,string information) const{
    if(!condition)
        throw(runtime_error("WadipolesT::Fail_"+information));
}
void WadipolesT::log(string information) const{
    GLOG->stream() << " " << information << endl;
}
