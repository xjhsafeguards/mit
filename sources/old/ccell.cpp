#include "ccell.h"

Inputfile& operator>>(Inputfile& inf, Cell& cel){
    if(inf.type() == "cif") cel.read_cif(inf.stream());
    else
        cel.test(0,"Cell::operator<<::wrong_inputfile_type:" + inf.type());
    return inf;
}
Outputfile& operator<<(Outputfile& of, const Cell& cel){
    if(of.type() == "in") cel.write_in(of.stream());
    else if(of.type() == "POSCAR") cel.write_POSCAR(of.stream());
    else
        cel.test(0,"Cell::operator<<::wrong_outputfile_type:" + of.type());
    return of;
}
Outputfile& operator<<(Outputfile& of,Cell& cel){
    if(of.type() == "in") cel.write_in(of.stream());
    else if(of.type() == "POSCAR") cel.write_POSCAR(of.stream());
    else
        cel.test(0,"Cell::operator<<::wrong_outputfile_type:" + of.type());
    return of;
}

unsigned Cell::start=1;
unsigned Cell::step=1;
unsigned Cell::end=-1;
map<string,int> Cell::index_dic;

void Cell::init(){
    os() << "Read in Cells" << endl;
    for(int i=1;i<start;i++)
        skip();
    set();
}
bool Cell::next(){
    if(num == 0 )
        init();
    else if(!inrange()){
        return false;
    }else{
        for(int i=1;i<step;i++){
            skip();
        }
        set();
    }
    return true;
}

Cell::Cell(Input input): cellp(make_shared<Cellp>(input.cellp)),atoms(make_shared<vector<Atom> >(input.ntype)),wannier(make_shared<Wannier>(input.nband)),geo_file(input.geo_file),cel_file(input.cel_file),wan_file(input.wan_file){
    if(wan_file.allocate())
        test(input.nband>0,"Read_wannier::need_nband");
    if(geo_file.allocate()){
        test(input.ntype>0,"Read_position::need_ntype");
        test(input.ntype==input.atom_num.size(),"Read_position::need_correct_atom_num");
    }
    auto name_it=input.atom_name.cbegin();
    auto num_it=input.atom_num.cbegin();
    auto charge_it=input.atom_charge.cbegin();
    auto mass_it=input.atom_mass.cbegin();
    for(auto it=atoms->begin(); it!=atoms->end() and name_it!= input.atom_name.cend() ;it++,name_it++)
        it->set_name(*name_it);
    for(auto it=atoms->begin(); it!=atoms->end() and num_it!= input.atom_num.cend() ;it++,num_it++)
        it->set_na(*num_it);
    for(auto it=atoms->begin(); it!=atoms->end() and charge_it!= input.atom_charge.cend() ;it++,charge_it++)
        it->set_charge(*charge_it);
    for(auto it=atoms->begin(); it!=atoms->end() and mass_it!= input.atom_mass.cend() ;it++,mass_it++)
        it->set_mass(*mass_it);
    start = input.ss_start;
    step = input.ss_step;
    end = input.ss_stop;
}
int Cell::anum() const{
    int sum=0;
    for_each(atoms->cbegin(),atoms->cend(),[&sum](const Atom& at){sum += at.num();});
    return sum;
}
int Cell::anum(int i) const{
#ifndef N_TEST
    test(allocate_atoms(i),"Cell::anum::out_of_range");
#endif
    return atoms->at(i).num();
}
int Cell::anum(string name) const{
    return anum(aindex(name));
}
int Cell::ana() const{
    int sum=0;
    for_each(atoms->cbegin(),atoms->cend(),[&sum](const Atom& at){sum += at.get_na();});
    return sum;
}
int Cell::ana(int i) const{
#ifndef N_TEST
    test(allocate_atoms(i),"Cell::ana::out_of_range");
#endif
    return atoms->at(i).get_na();
}
int Cell::ana(string name) const{
    return ana(aindex(name));
}
int Cell::aindex(string name) const{
    if(index_dic.count(name) == 0){
        int i=0;
        for(auto it=atoms->cbegin();it!=atoms->cend();it ++)
        {
            if(it->get_name() == name)
                return i;
            else i++;
        }
        test(0, "Get_atom_index::No_match_name");
        index_dic.insert({name,i});
    }
    return index_dic[name];
}
int Cell::wnum() const{
#ifndef N_TEST
    test(allocate_wannier(),"Cell::wnum::out_of_range");
#endif
    return wannier->num();
}
int Cell::wna() const{
#ifndef N_TEST
    test(allocate_wannier(),"Cell::wnum::out_of_range");
#endif
    return wannier->get_na();
}
double Cell::avolume() const{
    return (count!=0 ? volume_sum/count : volume());
}
shared_ptr<Cell> Cell::save() const{
    shared_ptr<Cell> tmpcell = make_shared<Cell> (*this);
    tmpcell->cellp = cellp->save();
    tmpcell->atoms = make_shared<vector<Atom> >();
    for( const auto& at : *atoms){
        tmpcell->atoms->push_back(*(at.save()));
    }
    tmpcell->wannier = wannier->save();
    //tmpcell->fatoms = make_shared<vector<Atom> >();
    //for( const auto& at : *fatoms){
    //    tmpcell->fatoms->push_back(*(at.save()));
    //}
    //tmpcell->fwannier = fwannier->save();
    return tmpcell;
}
typename Cell::position_type Cell::Cartesian(const cellp_value_type& cp,const position_type& a){
    return a*cp;
}
typename Cell::position_type Cell::In_cell(const position_type& a){
    position_type b(a);
    b.normal_BC(position_type(1.0,1.0,1.0));
    return b;
}
typename Cell::position_type Cell::Shortest_vector(const position_type& a,const position_type& b){
    position_type c = b - a;
    c = In_cell(c+position_type(0.5,0.5,0.5)) - position_type(0.5,0.5,0.5);
    return c;
}
double Cell::Distance(const cellp_value_type& cp,const position_type& a,const position_type& b){
    return Cartesian(cp,Shortest_vector(a,b)).norm();
}
double Cell::distance(int i1,int j1, int i2, int j2) const{
    return cartesian(shortest_vector(afpos(i1,j1),afpos(i2,j2))).norm();
}
double Cell::distance(string n1,int j1,string n2, int j2) const{
    return distance(aindex(n1),j1,aindex(n2),j2);
}
double Cell::distance(int i1,int j1, int w2) const{
    return cartesian(shortest_vector(afpos(i1,j1),wfpos(w2))).norm();
}
double Cell::distance(string n1,int j1, int w2) const{
    return cartesian(shortest_vector(afpos(n1,j1),wfpos(w2))).norm();
}
double Cell::angle(int i1,int j1, int i2, int j2, int i3, int j3) const{
    vector_type v1=shortest_vector(afpos(i1,j1),afpos(i2,j2));
    vector_type v2=shortest_vector(afpos(i3,j3),afpos(i2,j2));
    v1 = cartesian(v1);
    v2 = cartesian(v2);
    return acos(v1*v2/v1.norm()/v2.norm())/PI*180;
}
double Cell::angle(string n1,int j1,string n2, int j2,string n3, int j3) const{
    return angle(aindex(n1),j1,aindex(n2),j2,aindex(n3),j3);
}
typename Cell::vector_type Cell::cvector(int i1,int j1, int i2, int j2) const{
    return cartesian(shortest_vector(afpos(i1,j1),afpos(i2,j2)));
}
typename Cell::vector_type Cell::cvector(string n1,int j1,string n2, int j2) const{
    return cartesian(shortest_vector(afpos(n1,j1),afpos(n2,j2)));
}
typename Cell::vector_type Cell::cvector(int i1,int j1, int w2) const{
    return cartesian(shortest_vector(afpos(i1,j1),wfpos(w2)));
}
typename Cell::vector_type Cell::cvector(string n1,int j1, int w2) const{
    return cartesian(shortest_vector(afpos(n1,j1),wfpos(w2)));
}
void Cell::read_cif(istream& is){
    set_read_cif(is);
    cellp->read_cif(is);
    for( auto& at: (*atoms)){
        at.read_cif(is);
    }
    set_up_fatoms();
}
void Cell::write_in(ostream& os) const{
#ifndef N_TEST
    test(allocate_cellp(),"Cell::write_POSCAR::no_cellp");
    test(allocate_atoms(),"Cell::write_POSCAR::no_atoms");
    test(allocate_fatoms(),"Cell::write_POSCAR::no_atoms");
#endif
    os << "CELL_PARAMETERS {angstrom}" << endl;
    cellp->write_POSCAR(os);
    os << "ATOMIC_POSITIONS {crystal}" << endl;
    for( const auto& at: *fatoms)
        at.write_in(os);
}
void Cell::write_POSCAR(ostream& os) const{
#ifndef N_TEST
    test(allocate_cellp(),"Cell::write_POSCAR::no_cellp");
    test(allocate_atoms(),"Cell::write_POSCAR::no_atoms");
    test(allocate_fatoms(),"Cell::write_POSCAR::no_atoms");
#endif
    os << endl;
    os << 1.0 << endl;
    cellp->write_POSCAR(os);
    for( const auto& at: *atoms)
        os << " "<< at.get_name();
    os << endl << " ";
    for( const auto& at: *atoms)
        os << " "<< at.get_na();
    os << endl;
    os << "Direct" << endl;
    for( const auto& at: *fatoms)
        at.write_POSCAR(os);
}
void Cell::set_atoms_nfractional(){
    for( auto& atom: (*atoms) )
        set_pos_cartesian(atom);
}
void Cell::set_wannier_nfractional(){
    set_pos_cartesian(*wannier);
}
void Cell::set_fatoms_fractional(){
    for( auto& atom: (*fatoms) )
        set_pos_friction(atom);
}
void Cell::set_fwannier_fractional(){
    set_pos_friction(*fwannier);
}
void Cell::set_up_fatoms(){
    fatoms = value_copy(atoms);
    for(auto it=fatoms->begin();it!=fatoms->end();it++)
        it->value_duplicate();
    set_atoms_nfractional();
    set_fatoms_fractional();
}
void Cell::set_geo(){
    if(geo_file.allocate()){
        if(geo_file.type() == "pos"){
            (atoms->begin())->set_label_file(geo_file);
            check_labels(atoms->begin()->get_labels());
        }else
            test(0,"set_geo::wrong_geo_file_type");
        for(auto it=atoms->begin();it!=atoms->end();it++)
            geo_file >> *it;
        fatoms = value_copy(atoms);
        for(auto it=fatoms->begin();it!=fatoms->end();it++)
            it->value_duplicate();
        set_atoms_nfractional();
        set_fatoms_fractional();
    }
}
void Cell::set_wan(){
    if(wan_file.allocate()){
        wan_file >> *wannier;
        check_labels(wannier->get_labels());
        fwannier = value_copy(wannier);
        fwannier->value_duplicate();
        set_wannier_nfractional();
        set_fwannier_fractional();
    }
}
void Cell::set_cel(){
    if(cel_file.allocate()){
        cel_file >> *cellp;
        check_labels(cellp->get_labels());
    }
}
void Cell::skip_geo(){
    if(geo_file.allocate()){
        (atoms->begin())->set_label_file(geo_file);
        for(auto it=atoms->begin();it!=atoms->end();it++)
            it->skip(geo_file);
    }
}
void Cell::skip_wan(){
    if(wan_file.allocate())
        wannier->skip(wan_file);
}
void Cell::skip_cel(){
    if(cel_file.allocate())
        cellp->skip(cel_file);
}
void Cell::set_pos_friction(Position_data& pos){
    if(!pos.isfractional()){
        cellp_value_type cp(box().inverse());
        for(auto it=pos.begin();it!=pos.end();it++){
            *it = (*it)*cp;
        }
        pos.set_fractional(true);
    }
}
void Cell::set_pos_cartesian(Position_data& pos){
    if(pos.isfractional()){
        cellp_value_type cp(box());
        for(auto it=pos.begin();it!=pos.end();it++){
            *it = (*it)*cp;
        }
        pos.set_fractional(false);
    }
}
void Cell::check_labels(vector<string> label){
    if(label.size() == 0)
        return;
    if(!allocate_label())
        set_label(label);
    else
        for(int i=0; i<2; i++){
            if((*labels)[i] != label[i])
                cout << "Fail: "<< (*labels)[i] << " and " << label[i] << endl;
            test((*labels)[i] == label[i],"Cell::check_labels:" + to_string(i));
        }
}
void Cell::set_read_cif(istream& is){
    is.seekg(0,is.beg);
    cellp = make_shared<Cellp>();
    atoms = make_shared<vector<Atom>>();
    string buffer;
    Seek_value(is,"_chemical_formula_sum");
    getline(is,buffer);
    string name,na;
    for(auto letter : buffer){
        if(isalpha(letter))
            name.push_back(letter);
        else if(isdigit(letter))
            na.push_back(letter);
        else if( name.size() and na.size() ){
            log(name + ":" + na);
            atoms->push_back(Atom(name,stoi(na)));
            name.clear();
            na.clear();
        }
    }
}
/*
 Cell::Cell(int ntype):atoms(ntype){
 }
 Cell::Cell(int ntype, int nwan): atoms(ntype), wannier(nwan) {
 }
 */
