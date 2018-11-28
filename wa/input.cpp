#include "input.h"
#include "anwater.h"
#include "wahbond.h"
#include "wadipole.h"

Input INPUT;

Inputfile& operator>>(Inputfile& inf, Input& input){
    input.log("Read "+inf.name()+" for input");
    input.log() << endl;
    if(inf.type() == "INPUT") input.read_input(inf.stream());
    else
        input.test(0,"INPUT::operator>>::wrong_inputfile_type:" + inf.type());
    input.log() << endl;
    return inf;
}

void Input::read_input(istream& is){
    test_read(is);
    
    while(is >> buffer)
    {
        read = false;
        if(buffer[0]!='#' and buffer[0]!='!'){
            read_value("calculation",is,calculation);
            read_value("type",is,type);
            read_value("ensemble",is,ensemble);
            read_value("system",is,system);
            read_value("geo_file",is,geo_file);
            read_value("cel_file",is,cel_file);
            read_value("wan_file",is,wan_file);
            read_values("in_file",is,in_files);
            read_value("out_file_type",is,out_files_type);
            //read_values("out_file",is,out_files);
            read_value("ss_start",is,ss_start);
            read_value("ss_stop",is,ss_stop);
            read_value("ss_step",is,ss_step);
            read_value("ss_n",is,ss_n);
            read_uvalue("dt",is,dt);
            read_value("ntype",is,ntype);
            read_value("nband",is,nband);
            read_values("atom_num",is,atom_num);
            read_values("atom_name",is,atom_name);
            read_values("atom_mass",is,atom_mass);
            read_values("atom_charge",is,atom_charge);
            read_uvalue("cellp",is,cellp);
            read_uvalue("T",is,T);
            read_uvalue("OH_distance",is,Water::OH_distance);
            read_uvalue("OO_distance",is,Wahbonds::OO_distance);
            read_uvalue("HOO_angle",is,Wahbonds::HOO_angle);
            read_uvalue("OW_distance",is,Wadipoles::OW_distance);
            //read_uvalue("HW_distance",is,HW_distance);
            read_uvalue("cutoff",is,cutoff);
            read_uvalue("alpha",is,alpha);
            read_value("thicken",is,thicken);
            read_uvalue("upper_limit",is,upper_limit);
            read_value("delta",is,delta);
            read_uvalue("eps",is,eps);
            read_uvalue("displacement",is,displacement);
            read_value("tcf_max",is,tcf_max);
            read_values("index",is,index);
            read_values("parameter",is,parameter);
        }
        buffer = "#";
        if(!read)
            is.ignore(500,'\n');
    }
}

int Input::get_atom_num() const{
    return accumulate(atom_num.cbegin(),atom_num.cend(),0);
}
int Input::get_atom_num(size_t i) const{
#ifndef N_TEST
    test(i<atom_num.size(),"Get_atom_num::Out_of_range");
#endif
    return atom_num.at(i);
}
int Input::get_atom_num(string name) const{
    return get_atom_num(get_atom_index(name));
}
size_t Input::get_atom_index(string name) const{
    auto it = find(atom_name.cbegin(), atom_name.cend(), name);
    test(it!= atom_name.cend(), "Get_atom_index::No_match_name");
    return it - atom_name.cbegin();
}

void Input::read_values(string key, istream& is, vector<Inputfile>& data){
    if(buffer == key and !read){
        test_read(is);
        getline(is,buffer);
        istringstream iss(buffer);
        log() << setw(15) << key << endl;
        while(iss >> buffer)
            data.push_back(buffer);
        read = true;
    }
}
void Input::read_uvalue(string key, istream& is, Cellp& data){
    if(buffer == key and !read){
        test_read(is);
        getline(is,buffer);
        istringstream iss(buffer);
        iss >> unit;
        is >> data;
        data() = data()*Unit::unitconv(unit);
        unit = 1;
        log() << setw(15) << key << '\n';
        for(int i=0; i!= 3; i++)
            log() << setw(15) << data[i] << endl;
        read = true;
    }
}
