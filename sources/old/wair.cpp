#include "wair.h"
#include "input.h"

void WaIR::Routine(){
    Help::help("wair",INPUT.type);
    switch(INPUT.type){
        case 0:
            WaIR::C0();
            break;
        case 1:
            WaIR::C80();
            break;
        case 2:
            WaIR::C81();
            break;
        case 3:
            WaIR::C82();
            break;
        case 4:
            WaIR::C90();
            break;
        default:
            WaIR::TEST();
            break;
    }
}

void WaIR::C0(){
    Outputfile of("vdipole.txt");
    Outputfile of1("TCF.txt");
    Outputfile of2("IR.txt");
    Cell cel(INPUT);
    WadipolesT WDT;
    while(cel.next()){
        Wadipoles WD(cel);
        WD.set_dipoles_dipole();
        WDT.add(WD);
    }
    WDT.calculate_vdipole();
    WDT.write_vdipole(of.stream());
    IR ir;
    ir.set_dt(WDT.get_dt());
    ir.set_V(cel.avolume());
    ir.set_vdipole(WDT.get_vdipole());
    ir.calculate_ir();
    ir.write_tcf(of1.stream());
    ir.write_ir(of2.stream());
}

void WaIR::C80(){
    Outputfile of("vdipole.txt");
    Outputfile ofcp("cellp.txt");
    Outputfile of1("TCF.txt");
    Outputfile of2("IR.txt");
    Cell cel(INPUT);
    int nO = cel.ana("O");
    vector<typename Cell::cellp_value_type> cellpT;
    vector<vector<typename Cell::position_type> > dposT(nO);
    vector<WadipolesT> vWDT(nO);
    while(cel.next()){
        Wadipoles WD(cel);
        WD.set_dipoles_dipole();
        double Time = WD.time();
        cellpT.push_back(cel.box());
        for(int i=0;i<nO;i++){
            vWDT[i].add(WD.get_dipole(i),Time);
            dposT[i].push_back(WD.fposO(i));
        }
    }
    IR3<tuple<int,int>> ir;
    ir.set_dt(vWDT[0].get_dt());
    ir.set_V(cel.avolume());
    for( auto& cp : cellpT){
        ofcp.stream() << cp << endl;
    }
    for( auto& WDT : vWDT){
        WDT.calculate_vdipole();
        WDT.write_vdipole(of.stream());
        ir.push_back_vdipole(WDT.get_vdipole());
    }
    ir.label_condition();
    // get<0>(m) index m,get<0>(n),get<1>(m) snapshot m,get<1>(n)
    ir.calculate_tcf_all([&](tuple<int,int> m,tuple<int,int> n) -> bool {
        int im=get<0>(m),sm=get<1>(m)+1,in=get<0>(n),sn=get<1>(n)+1;
        return ( im!=in and WaIR::in_range(Cell::Distance(cellpT.at(sn),dposT.at(im).at(sm),dposT.at(in).at(sn))));
    });
    ir.calculate_ir();
    ir.write_tcf(of1.stream());
    ir.write_ir(of2.stream());
}

void WaIR::C81(){
    Outputfile of("vdipole.txt");
    Cell cel(INPUT);
    int nO = cel.ana("O");
    vector<Wadipoles> WDTS;
    vector<WadipolesT> vWDT(nO);
    vector<typename Cell::cellp_value_type> cellpT;
    vector<vector<typename Cell::position_type> > dposT(nO);
    while(cel.next()){
        Wadipoles WD(cel);
        WD.set_hbonds();
        WD.set_dipoles_dipole();
        double Time = WD.time();
        cellpT.push_back(cel.box());
        for(int i=0;i<nO;i++){
            vWDT[i].add(WD.get_dipole(i),Time);
            dposT[i].push_back(WD.fposO(i));
        }
        WD.unset_dipoles();
        WD.unset_waters();
        WD.unset_cel();
        WDTS.push_back(WD);
    }
    IR3<tuple<int,int>> ir;
    ir.set_dt(vWDT[0].get_dt());
    ir.set_V(cel.avolume());
    for( auto& WDT : vWDT){
        WDT.calculate_vdipole();
        WDT.write_vdipole(of.stream());
        ir.push_back_vdipole(WDT.get_vdipole());
    }
    ir.label_condition();
    for( auto& nu : INPUT.index){
        Outputfile of1("TCF_" + to_string(nu) +"_Hb.txt");
        Outputfile of2("IR_" + to_string(nu) +"_Hb.txt");
        // get<0>(m) index m,get<0>(n),get<1>(m) snapshot m,get<1>(n)
        ir.calculate_tcf_all([&](tuple<int,int> m,tuple<int,int> n) -> bool {
            int im=get<0>(m),sm=get<1>(m)+1,in=get<0>(n),sn=get<1>(n)+1;
            //if( im!=in and WaIR::in_range(Cell::Distance(WDTS[sn].get_cell().box(),WDTS[sm].fposO(im),WDTS[sn].fposO(in))) ){
            if( im!=in and WaIR::in_range(Cell::Distance(cellpT.at(sn),dposT.at(im).at(sm),dposT.at(in).at(sn))) ){
                Wadipoles& WD = WDTS[sn];
                int bonded = 0;
                for( auto& ind : WD.get_index_accept(im) ){
                    if( WD.check_bond(in,ind) )
                        ++bonded;
                }
                for( auto& ind : WD.get_index_donate(im) ){
                    if( WD.check_bond(in,ind) )
                        ++bonded;
                }
                if(bonded==nu){
                    return true;
                }
            }
               return false;
        });
        ir.calculate_ir();
        ir.write_tcf(of1.stream());
        ir.write_ir(of2.stream());
        ir.unset_tcf();
    }
}

void WaIR::C82(){
    Outputfile of("vdipole.txt");
    Cell cel(INPUT);
    int nO = cel.ana("O");
    vector<Wadipoles> WDTS;
    vector<WadipolesT> vWDT(nO);
    vector<typename Cell::cellp_value_type> cellpT;
    vector<vector<typename Cell::position_type> > dposT(nO);
    while(cel.next()){
        Wadipoles WD(cel);
        WD.set_hbonds();
        WD.set_dipoles_dipole();
        double Time = WD.time();
        cellpT.push_back(cel.box());
        for(int i=0;i<nO;i++){
            vWDT[i].add(WD.get_dipole(i),Time);
            dposT[i].push_back(WD.fposO(i));
        }
        WD.unset_dipoles();
        WD.unset_waters();
        WD.unset_cel();
        WDTS.push_back(WD);
    }
    IR3<tuple<int,int>> ir;
    ir.set_dt(vWDT[0].get_dt());
    ir.set_V(cel.avolume());
    for( auto& WDT : vWDT){
        WDT.calculate_vdipole();
        WDT.write_vdipole(of.stream());
        ir.push_back_vdipole(WDT.get_vdipole());
    }
    ir.label_condition();
    for( auto& nu : INPUT.index){
        Outputfile of1("TCF_" + to_string(nu) +"_Hbc.txt");
        Outputfile of2("IR_" + to_string(nu) +"_Hbc.txt");
        // get<0>(m) index m,get<0>(n),get<1>(m) snapshot m,get<1>(n)
        ir.calculate_tcf_all([&](tuple<int,int> m,tuple<int,int> n) -> bool {
            int im=get<0>(m),sm=get<1>(m)+1,in=get<0>(n),sn=get<1>(n)+1;
            if( im!=in and WaIR::in_range(Cell::Distance(cellpT.at(sn),dposT.at(im).at(sm),dposT.at(in).at(sn))) ){
                Wadipoles& WD = WDTS[sn];
                int bonded = 0;
                if( WD.check_bond(in,im) )
                    ++bonded;
                if(bonded==nu){
                    return true;
                }
            }
            return false;
        });
        ir.calculate_ir();
        ir.write_tcf(of1.stream());
        ir.write_ir(of2.stream());
        ir.unset_tcf();
    }
}

void WaIR::C90(){
    Cell cel(INPUT);
    int nO = cel.ana("O");
    vector<Wadipoles> WDTS;
    vector<WadipolesT> vWDT(nO);
    while(cel.next()){
        Wadipoles WD(cel);
        WD.set_hbonds();
        WD.set_dipoles_dipole();
        double Time = WD.time();
        for(int i=0;i<nO;i++){
            vWDT[i].add(WD.get_dipole(i),Time);
        }
        WD.unset_dipoles();
        WD.unset_waters();
        WD.unset_cel();
        WDTS.push_back(WD);
    }
    Outputfile of("vdipole.txt");
    IR3<tuple<int,int>> ir;
    ir.set_dt(vWDT[0].get_dt());
    ir.set_V(cel.avolume());
    for( auto& WDT : vWDT){
        WDT.calculate_vdipole();
        WDT.write_vdipole(of.stream());
        ir.push_back_vdipole(WDT.get_vdipole());
    }
    ir.label_condition();
    
    for( auto& nu : INPUT.index){
        Outputfile of1("TCF_" + to_string(nu) +"_Hb.txt");
        Outputfile of2("IR_" + to_string(nu) +"_Hb.txt");
        // get<0>(m) index m,get<0>(n),get<1>(m) snapshot m,get<1>(n)
        ir.calculate_tcf_all([&](tuple<int,int> m,tuple<int,int> n) -> bool {
            int im=get<0>(m),sm=get<1>(m)+1,in=get<0>(n),sn=get<1>(n)+1;
            if( im!=in and WDTS[sn].nhbond(im)==nu and WDTS[sn].check_bond(in,im) )
                    return true;
            return false;
        });
        ir.calculate_ir();
        ir.write_tcf(of1.stream());
        ir.write_ir(of2.stream());
        ir.unset_tcf();
    }
}

void WaIR::TEST(){
    Inputfile of("vdipole.txt");
    Outputfile of1("TCF.txt");
    Outputfile of2("IR.txt");
    Cell cel(INPUT);
    int nO = cel.ana("O");
    vector<WadipolesT> vWDT(nO);
    IR3<tuple<int,int>> ir;
    for( auto& WDT : vWDT){
        WDT.read_vdipole(of.stream());
        ir.push_back_vdipole(WDT.get_vdipole());
    }
    ir.set_dt(vWDT[0].get_dt());
    ir.set_V(cel.avolume());
    ir.label_condition();
    // get<0>(m) index m,get<0>(n),get<1>(m) snapshot m,get<1>(n)
    ir.calculate_tcf_all([&](tuple<int,int> m,tuple<int,int> n) -> bool {
        int im=get<0>(m),sm=get<1>(m),in=get<0>(n),sn=get<1>(n);
        return ( im == in );
    });
    ir.calculate_ir();
    ir.write_tcf(of1.stream());
    ir.write_ir(of2.stream());
}
