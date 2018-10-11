#include "ccellp.h"

Inputfile& operator>>(Inputfile &inf,Cellp &celp){
    if(inf.type() == "cif") celp.read_cif(inf.stream());
    else if(inf.type() == "cel") celp.read_cel(inf.stream());
    else
        celp.test(0,"CELLP::operator>>::wrong_inputfile_type:" + inf.type());
    return inf;
}
Outputfile& operator<<(Outputfile& of, const Cellp& celp){
    if(of.type() == "in") celp.write_in(of.stream());
    else if(of.type() == "POSCAR") celp.write_POSCAR(of.stream());
    else
        celp.test(0,"CELLP::operator<<::wrong_outputfile_type:" + of.type());
    return of;
}

Cellp::Cellp(value_type in_mat): mat(new value_type(in_mat)){
    
}
shared_ptr<Cellp> Cellp::save() const{
    shared_ptr<Cellp> tmpcp = make_shared<Cellp>(*this);
    tmpcp->mat = value_copy(mat);
    return tmpcp;
}
void Cellp::get_parameter(double& a,double& b,double& c,double& al,double& be,double& ga) const{
    a = (*mat)[0].norm();
    b = (*mat)[1].norm();
    c = (*mat)[2].norm();
    al = (*mat)[1].vangle((*mat)[2]);
    be = (*mat)[2].vangle((*mat)[0]);
    ga = (*mat)[0].vangle((*mat)[1]);
}
void Cellp::write_mat(ostream& os,double unitconv) const{
    auto buffer = *mat * unitconv;
    os << showpoint << fixed << buffer << noshowpoint << endl;
    os.unsetf(std::ios_base::floatfield);
}
void Cellp::set_mat(istream& is,double unitconv){
    value_type buffer;
    is >> buffer;
    mat = make_shared<value_type>(buffer*unitconv);
}
void Cellp::skip(Inputfile& inf) const{
    if(inf.type() == "cel") skip_cel(inf.stream());
    else
        test(0,"CELLP::skip::wrong_inputfile_type:" + inf.type());
}
void Cellp::set_parameter(double a,double b,double c,double al,double be,double ga,double unitconv){
    mat = make_shared<value_type>();
    a=a*unitconv; b=b*unitconv; c=c*unitconv;
    al = al * PI /180;
    be = be * PI /180;
    ga = ga * PI /180;
    double cosga = cos(ga);
    double singa = sin(ga);
    double cosbe = cos(be);
    double cosal = cos(al);
    (*mat)[0].set(a,0,0);
    (*mat)[1].set(b*cosga,b*singa,0);
    (*mat)[2].set(c*cosbe,c*(cosal-cosbe*cosga)/singa,c*sqrt(1+2*cosal*cosbe*cosga-cosal*cosal-cosbe*cosbe-cosga*cosga)/singa);
}
void Cellp::read_cif(istream& is){
    is.seekg(0,is.beg);
    double a,b,c,al,be,ga;
    Seek_value(is,"_cell_length_a",a);
    Seek_value(is,"_cell_length_b",b);
    Seek_value(is,"_cell_length_c",c);
    Seek_value(is,"_cell_angle_alpha",al);
    Seek_value(is,"_cell_angle_beta",be);
    Seek_value(is,"_cell_angle_gamma",ga);
    set_parameter(a,b,c,al,be,ga);
}

istream& operator>>(istream& is,Cellp& cp){
    typename Cellp::value_type buffer;
    is >> buffer;
    cp.mat = make_shared<Cellp::value_type>(buffer);
    return is;
}
ostream& operator<<(ostream& os,const Cellp& cp){
    cp.test(cp.allocate_mat(),"Cellp::operator<<::no_mat");
    return os << cp.mat;
}

