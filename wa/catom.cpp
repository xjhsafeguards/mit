#include "catom.h"

Atom::Atom(string in_name, int in_na, double in_mass, double in_charge): Position_data(in_name,in_na), mass(in_mass), charge(in_charge){
}

Inputfile& operator>>(Inputfile &inf,Atom &at){
    if(inf.type() == "pos") at.read_pos(inf.stream());
    else if(inf.type() == "cif") at.read_cif(inf.stream());
    else
        at.test(0,"ATOM::operator>>::wrong_inputfile_type:" + inf.type());
    return inf;
}
Outputfile& operator<<(Outputfile &of, const Atom &at){
    if(of.type() == "in") at.write_in(of.stream());
    else if(of.type() == "POSCAR") at.write_POSCAR(of.stream());
    else
            at.test(0,"ATOM::operator<<::wrong_outputfile_type:" + of.type());
    return of;
}
shared_ptr<Atom> Atom::save() const{
    shared_ptr<Atom> tmpatom = make_shared<Atom>(*this);
    tmpatom->labels = value_copy(labels);
    tmpatom->position = value_copy(position);
    tmpatom->index = value_copy(index);
    return tmpatom;
}
void Atom::write_in(ostream& os) const{
    for( const auto& pos : *position)
        os << name << " " << setw(15) << pos[0] << " " << setw(15) << pos[1] << " " << setw(15) << pos[2] << endl;
}
void Atom::write_in(ostream& os, double unitconv) const{
    for( const auto& pos : *position){
        os << name << " " << setw(15) << pos[0] * unitconv  << " " << setw(15) << pos[1] * unitconv << " " << setw(15) << pos[2] * unitconv  << endl;
    }
}
void Atom::write_POSCAR(ostream& os) const{
    for( const auto& pos : *position)
        os << " " << setw(15) << pos[0]  << " " << setw(15) << pos[1]  << " " << setw(15) << pos[2] << endl;
}
void Atom::write_POSCAR(ostream& os, double unitconv) const{
    for( const auto& pos : *position)
        os << " " << setw(15) << pos[0] * unitconv  << " " << setw(15) << pos[1] * unitconv << " " << setw(15) << pos[2] * unitconv << endl;
}
void Atom::skip(Inputfile& inf) const{
    if(inf.type() == "pos") skip_pos(inf.stream());
    else
        test(0,"ATOM::skip::wrong_inputfile_type:" + inf.type());
}
void Atom::set_label_file(Inputfile& inf){
    if(inf.type() == "pos") set_label(inf.stream());
    else
        test(0,"ATOM::set_label::wrong_inputfile_type:" + inf.type());
}

void Atom::read_pos(istream& is){
    try{
        fractional = false;
        set_position(is,Unit::Bohr2A);
    }catch(runtime_error err){
        further_throw(err,"::ATOM:" + name + "::read_pos");
    }
}
void Atom::skip_pos(istream& is) const{
    try{
        skip_position(is);
    }catch(runtime_error err){
        further_throw(err,"::ATOM:" + name + "::skip_pos");
    }
}
void Atom::read_cif(istream& is){
    
    is.seekg(0,is.beg);
    string buffer;
    bool read_lable = false;
    bool fract = true;
    int count = 0;
    vector<string> atom_site;
    while(is){
        Seek_value(is,"loop_");
        is >> buffer;
        //cout << buffer;
        if(buffer.find("_atom_site_")!=string::npos){
            do{
                atom_site.push_back(buffer);
                is >> buffer;
                //cout << buffer;
            }while(buffer.find("_atom_site_")!=string::npos);
            position_type temp;
            while(count < get_na()){
                for( const auto& asite : atom_site){
                    if(asite == "_atom_site_type_symbol"){
                        if(buffer == get_name())
                            read_lable = true;
                        else
                            read_lable = false;
                    }
                    if(asite == "_atom_site_fract_x" and read_lable){
                        temp[0] = stod(buffer);
                    }
                    if(asite == "_atom_site_fract_y" and read_lable){
                        temp[1] = stod(buffer);
                    }
                    if(asite == "_atom_site_fract_z" and read_lable){
                        temp[2] = stod(buffer);
                        position->push_back(temp);
                        count ++;
                    }
                    if(is){
                        is >> buffer;
                    }
                    else
                        throw(runtime_error("Read_cif::can't_find_enough_"+ get_name()));
                }
            }
            set_fractional(fract);
            return;
        }
    }
    throw(runtime_error("Read_cif::no_atom_site"));
}

