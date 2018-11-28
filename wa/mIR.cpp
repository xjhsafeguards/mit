#include "mIR.h"
#include "input.h"

IR::IR(): T(INPUT.T),screen_t(INPUT.cutoff),screen_alpha(INPUT.alpha),tcf_nmax(INPUT.delta),ft_kmax(INPUT.upper_limit){
}
double IR::get_dt() const{
    test(allocate_dt(),"get_dt::no_dt");
    return dt;
}
const vector<double>& IR::get_tcf() const{
    test(allocate_tcf(),"get_tcf::no_tcf");
    return tcf;
}
void IR::write_tcf(ostream& os) const{
    int i = 0;
    for(const auto& tc : tcf)
        os << setw(20) << i++ * dt << " " << setw(20) << tc << endl;
}
void IR::write_ir(ostream& os) const{
    int i = 0;
    for(const auto& f : ir)
        os << setw(20) << i++ * dw << " " << setw(20) << f << endl;
}
void IR::unset_tcf(){
    tcf.clear();
    unset_ir();
}
void IR::unset_ir(){
    ir.clear();
}
void IR::set_screen(double st,double sa){
    screen_t = st;
    screen_alpha = sa;
}
void IR::set_dt(double in_dt){
    dt = in_dt;
}
void IR::set_T(double st){
    T = st;
}
void IR::set_V(double sv){
    V = sv;
}
void IR::set_dipole(shared_ptr<vector<dipole_type> > in_d){
    dipole = in_d;
}
void IR::set_vdipole(shared_ptr<vector<dipole_type> > in_v){
    vdipole = in_v;
}
void IR::set_tcf(const vector<double>& in_tcf){
    tcf = in_tcf;
}
void IR::set_up(const IR& in_ir){
    dt = in_ir.get_dt();
    tcf = in_ir.get_tcf();
}
void IR::calculate_vdipole(){
    test(allocate_dt(),"calculate_vdipole::no_dt"); 
    test(allocate_dipole(),"calculate_vdipole::no_dipole");
    c_vdipole(dipole,vdipole);
}
void IR::calculate_tcf(){
    if(!allocate_vdipole())
        calculate_vdipole();
    c_tcf_init(vdipole,vdipole,tcf_nmax);
    //correlate_function(vdipole->cbegin(),vdipole->cend(),back_inserter(tcf),nmax,1);
}
void IR::calculate_ir(){
    if(!allocate_tcf())
        calculate_tcf();
    log("Calculate ir with screen: t_cutoff:" + to_string(screen_t) + ", alpha:" + to_string(screen_alpha));
    gaussian_screen_tcf();
    //cout << setw(10) << dt << setw(10) <<kmax << setw(10) << tcf.size() << endl;
    int Ndx = 2*PI*tcf.size()/(dt*ft_kmax); // N == kmax/2PI*(dx*Ndx)
    dw = fourier_transform_function(tcf.cbegin(),tcf.cend(),back_inserter(ir),dt,Ndx,ft_kmax);
    unit_convert_ir();
}
void IR::c_vdipole(shared_ptr<vector<dipole_type> > d, shared_ptr<vector<dipole_type> > v){
    v = make_shared<vector<dipole_type> >();
    derivate(d->cbegin(),d->cend(),back_inserter(*v),dt);
}
void IR::c_tcf_init(shared_ptr<vector<dipole_type> > v1, shared_ptr<vector<dipole_type> > v2, int nmax){
    tcf.clear();
    correlate2_function(v1->cbegin(),v1->cend(),v2->cbegin(),v2->cend(),back_inserter(tcf),nmax,1);
}
void IR::c_tcf(shared_ptr<vector<dipole_type> > v1, shared_ptr<vector<dipole_type> > v2, int nmax){
    vector<double> tmp_tcf;
    correlate2_function(v1->cbegin(),v1->cend(),v2->cbegin(),v2->cend(),back_inserter(tmp_tcf),nmax,1);
    auto tmp = tcf.begin();
    for( auto i : tmp_tcf)
        *(tmp++) += i;
}
void IR::gaussian_screen_tcf(){
    test(allocate_dt(),"gaussian_screen_tcf::no_dt");
    test(allocate_tcf(),"gaussian_screen_tcf::no_tcf");
    double t = 0;
    double tmax = dt * (tcf.size() - 1);
    auto it = tcf.begin();
    for(; t<=screen_t; t+=dt)
        it++;
    for(; t<=tmax; t+=dt)
        *(it++) *= exp(-0.5*pow(screen_alpha,2)*pow((t-screen_t)/(tmax-screen_t),2));
}
void IR::unit_convert_ir(){
    test(T>0,"unit_convert_ir::wrong_T");
    test(V>0,"unit_convert_ir::wrong_V");
    dw /= Unit::Cm_12Ps_1; // from 2PI*ps-1 to cm-1
    double uc = Unit::Miu*Unit::C/(3*Unit::Kb*T*V*Unit::A3)*Unit::EA*Unit::EA/Unit::Ps*Unit::M2C;
    for(auto& f: ir)
        f *= uc;
}
void IR::test(bool condition,string information) const{
    if(!condition)
        throw(runtime_error("IR::Fail_"+information));
}
//void IR::log(string information) const{
//    ofs_log << " " << information << endl;
//}

//IR2
IR2::IR2(){
}
void IR2::set_dipole(const vector<shared_ptr<vector<dipole_type> > > &in_d){
    dipole = in_d;
}
void IR2::push_back_dipole(shared_ptr<vector<dipole_type> >  in_d){
    dipole.push_back(in_d);
}
void IR2::push_back_dipole(const vector<dipole_type>&  in_d){
    dipole.push_back(make_shared<vector<dipole_type> >(in_d));
}
void IR2::set_vdipole(const vector<shared_ptr<vector<dipole_type> > > &in_v){
    vdipole = in_v;
}
void IR2::push_back_vdipole(shared_ptr<vector<dipole_type> >  in_v){
    vdipole.push_back(in_v);
}
void IR2::push_back_vdipole(const vector<dipole_type>&  in_v){
    vdipole.push_back(make_shared<vector<dipole_type> >(in_v));
}
void IR2::calculate_vdipole(){
    test(allocate_dt(),"calculate_vdipole::no_dt");
    test(allocate_dipole(),"calculate_vdipole::no_dipole");
    vdipole.clear();
    for( auto &d : dipole){
        shared_ptr<vector<dipole_type> > v;
        c_vdipole(d,v);
        vdipole.push_back(v);
    }
}
void IR2::calculate_tcf(int i, int j){
    if(!allocate_vdipole(i) or !allocate_vdipole(j))
        calculate_vdipole();
    test(allocate_vdipole(i) and allocate_vdipole(j),"calculate_tcf::no_vdipole_i_j");
    if(allocate_tcf())
        c_tcf(vdipole.at(i),vdipole.at(j),tcf_nmax);
    else
        c_tcf_init(vdipole.at(i),vdipole.at(j),tcf_nmax);
}
//IR3
template<>
void IR3<tuple<int,int>>::label_condition(){
    test(allocate_vdipole(),"label_condition::no_vdipole");
    for(int j=0; j<vdipole.size(); j++){
        shared_ptr<vector<tuple<int,int>> > tmpc = make_shared<vector<tuple<int,int>> >();
        for(int i=0; i<vdipole[0]->size(); i++)
            tmpc->push_back(make_tuple(j,i));
        condition.push_back(tmpc);
    }

}
