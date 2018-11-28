#include "gdata.h"

void Data::set_label(istream &is){
    labels->clear();
    string s;
    label_type buffer;
    getline(is, s);
    istringstream ss(s);
    while(ss >> buffer){
        labels->push_back(buffer);
    }
}
void Data::set_label(vector<label_type> in_label){
    labels = make_shared<vector<label_type>>(in_label);
}
void Data::set_label_file(Inputfile &inf){
    set_label(inf.stream());
}

void Data::check_label(istream &is){
    if(allocate_label()){
        stringstream ss;
        string s;
        label_type buffer;
        getline(is, s);
        int i=0;
        while(ss){
            ss >> buffer;
            test(buffer==labels->at(i++),"Check_label");
        }
    }
    else
        set_label(is);
}

// Position_data
Position_data::Position_data(string in_name, int in_na): name(in_name),na(in_na) {
}

void  Position_data::set_position(istream &is,double unit_conv){
    unset_position();
    for(int i=0; i<get_na(); i++){
        test(is.good(),"set_position::istream_bad");
        position_push_back(is,unit_conv);
    }
}

void Position_data::value_duplicate(){
    position = value_copy(position);
}


//Dictionary_data



