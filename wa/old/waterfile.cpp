#include "waterfile.h"
#include "input.h"

Waterfile::Waterfile()
{
    O = -1;
    H = -1;
    allocate_waters = false;
    nwater =0;
    nion =0;
}

Waterfile::~Waterfile()
{
    if(allocate_waters)
        delete [] waters;
}

void Waterfile::Read_water(Cellfile &Cel)
{
    assert(Cel.allocate_atoms); //make sure Cel is not empty
    
    for(int i=0; i<Cel.ntype;i++) //find out index of O and H in Cel
    {
        if(Cel.atoms[i].id == "O") O=i;
        if(Cel.atoms[i].id == "H") H=i;
        if(Cel.atoms[i].id == "C") C=i;
    }
    if( O == -1 or H == -1)
    {
        cout << "Can not find O or H to Read_water" << endl;
        exit(1);
    }
        
    allocate_waters = true; //allocate positions for waters
    cel = &Cel; // mark the cellfile
    nwater = Cel.atoms[O].na;
    waters = new Water[nwater];
    
    for(int i=0; i<nwater ; i++) //start to find out waters
    {
        waters[i].indexO = i;
        waters[i].numH = 0;
        for(int j=0; j< Cel.atoms[H].na; j++)
        {
            double tmpdis;
            tmpdis = Cel.atoms[O].pos[i].distance_BC(Cel.atoms[H].pos[j],Cel.celldm); // considering boundary conditions
            if(tmpdis < INPUT.OH_distance/INPUT.unitconv)
            {
                //cout << i << '\t' << j << endl;
                waters[i].indexH[waters[i].numH] = j;
                waters[i].disH[waters[i].numH] = tmpdis;
                waters[i].numH++;
            }
        }
        //mark index of ions
        if(waters[i].numH != 2)
        {
            nion++;
            ions.push_back(i);
        }
    }
}


void Waterfile::Print_water(string out_file)
{
    cout << "Printing water of the System using type: " << INPUT.type << endl;

    ifstream ifs(INPUT.geo_file.c_str());
    ofstream ofs(out_file.c_str());
    ofstream ofs2("water_pos.txt");
    ofstream ofs3("ion_pos.txt");
    if(ifs.fail())
    {
        cerr << "Can not open in_file:" << INPUT.geo_file << endl;
        exit(1);
    }

    if(ofs3.fail())
    {
        cerr << "Can not open out_file: ion_count.txt to write!" << endl;
        exit(1);
    }
    if(ofs.fail())
    {
        cerr << "Can not open out_file:" << out_file << " to write!" << endl;
        exit(1);
    }
    if(ofs2.fail())
    {
        cerr << "Can not open out_file: water_pos.txt to write!" << endl;
        exit(1);
    }
    if(INPUT.type == 0 or INPUT.type == 1)
    {
        ofs << "index of O\t index of H\t index of H ...." << endl;
        ofs2 << "num of H\t index of O\t position of O\t" << endl;
    }
    if(INPUT.type == 0 or INPUT.type == 2)
    {
        ofs3 << "snapshot" << '\t' << "time" << '\t' << "index of ion" << '\t' << "numH in ion" << '\t' << "position x,y,z" << endl;
    }
    
    //Cellfile cel;
    for(int i=1; i<INPUT.ss_stop; i++)
    {
        Cellfile cel;
        
        if(i<INPUT.ss_start || (i-INPUT.ss_start)%INPUT.ss_step)
            cel.ReadGeometry(ifs,1); // skip the one we dont need
        else
        {
            if(cel.ReadGeometry(ifs)==0) break; // if to the end of the file skip
            cel.organize_pos();
            
            Waterfile WF;
            int count[8] = {0}; // use to count num of each kind
    
            WF.Read_water(cel);
            if(INPUT.type == 0 or INPUT.type == 1)
            {
                ofs << "snapshot :" << cel.snapshot << endl;
                ofs2 << "snapshot :" << cel.snapshot << endl;
                
                for(int i=0; i<WF.nwater ; i++)
                {
                    //writing ofs
                    ofs << WF.waters[i].indexO;
                    for(int j=0; j<WF.waters[i].numH; j++)
                        ofs << "    " << WF.waters[i].indexH[j];
                    ofs << endl;
                    count[WF.waters[i].numH]++;
                    //writing ofs2
                    ofs2 << WF.waters[i].numH <<"    " <<WF.waters[i].indexO << "    " <<cel.atoms[WF.O].pos[i].x*INPUT.unitconv << "    " <<cel.atoms[WF.O].pos[i].y*INPUT.unitconv << "    " << cel.atoms[WF.O].pos[i].z*INPUT.unitconv << endl;
                }
                for(int j=0; j<8; j++)
                    ofs << count[j] << "    ";
                ofs << endl;
            }
            if(INPUT.type == 0 or INPUT.type == 2)
            {
                 for(int j=0; j<WF.nion ; j++)
                 {   //writing of3
                     if(WF.waters[WF.ions[j]].numH !=2)
                         ofs3 << cel.snapshot << '\t' <<  cel.time << '\t' << WF.waters[WF.ions[j]].indexO  << '\t' << WF.waters[WF.ions[j]].numH << '\t' << cel.atoms[WF.O].pos[WF.ions[j]].x*INPUT.unitconv << "    " <<cel.atoms[WF.O].pos[WF.ions[j]].y*INPUT.unitconv << "    " << cel.atoms[WF.O].pos[WF.ions[j]].z*INPUT.unitconv << endl;
                 }
            }
        }
    }
    
}


void Waterfile::Add_constrain()
{
    cout << "Adding constrain on water to the System" << endl;
    
    ifstream ifs(INPUT.in_file.c_str());
    ofstream ofs("constrain.txt");
    if(ifs.fail())
    {
        cerr << "Can not open file: " << INPUT.in_file << " to read!" << endl;
        exit(1);
    }
    if(ofs.fail())
    {
        cerr << "Can not open file: constrain.txt to write" << endl;
        exit(1);
    }
    
    ofs << setprecision(8);
    string label;
    while(label != "ATOMIC_POSITIONS")
        Read_value(ifs,label); // find the position of atom positions
    
    Cellfile cel;
    Waterfile WF;
    
    cel.ReadInfile(ifs);
    WF.Read_water(cel);
    
    int nconstr = 0;
    for(int i=0; i<WF.nwater; i++)
    {
        if(WF.waters[i].numH == 1)
            nconstr += 1;
        if(WF.waters[i].numH == 2)
            nconstr += 3;
    }
    
    ofs << "CONSTRAINTS" << endl;
    ofs << nconstr << endl;
    
    
    
    for(int i=0; i<WF.nwater; i++)
    {
        if(WF.waters[i].numH == 1)
            ofs << "'distance'" << "      " << cel.atoms[WF.O].index[WF.waters[i].indexO] << "  " << cel.atoms[WF.H].index[WF.waters[i].indexH[0]] << endl;
        if(WF.waters[i].numH == 2)
        {
            ofs << "'distance'" << "      " << cel.atoms[WF.O].index[WF.waters[i].indexO] << "  " << cel.atoms[WF.H].index[WF.waters[i].indexH[0]] << endl;
            ofs << "'distance'" << "      " << cel.atoms[WF.O].index[WF.waters[i].indexO] << "  " << cel.atoms[WF.H].index[WF.waters[i].indexH[1]] << endl;
            ofs << "'planar_angle'" << "  " << cel.atoms[WF.H].index[WF.waters[i].indexH[0]] << "  " << cel.atoms[WF.O].index[WF.waters[i].indexO] << "  " << cel.atoms[WF.H].index[WF.waters[i].indexH[1]] << endl;
        }
    }
    
}

// add at 2018.2.2
Vector3<double> Waterfile::posO(int nO) const
{
    assert(allocate_waters);
    return cel->atoms[O].pos[waters[nO].indexO];
}

Vector3<double> Waterfile::posH(int nO, int nH) const
{
    my_throw(allocate_waters,"WANNIERFILE_POSWC: !allocate_waters");
    if(nH<waters[nO].numH)
        return cel->atoms[H].pos[waters[nO].indexH[nH]];
    else
    {
        throw(std::range_error("WATERFILE_POSH: nH >= numH"));
        return Vector3<double> (0,0,0);
        cerr << "Error, Can not read the " << nO << "th water's " << nH << "H position" << endl;
    }
}
// add at 2018.2.25
// updated for case 2 at 2018.3.16
void Waterfile::Animate_water()
{
    cout << "Animate water with type: " << INPUT.type << endl;
    switch (INPUT.type) {
        case 0:
            cout << "Making water vibration eigen mode animation" << endl;
            break;
        case 2:
            cout << "Bending water with displacement: " << INPUT.displacement << ", step: " << INPUT.delta << endl;
            break;
        default:
            cout << "Wrong INPUT type" << endl;
            break;
    }
    
    switch (INPUT.type){
        case 0:
        {
            ifstream ifs(INPUT.in_file.c_str());\
            file_assert(ifs,INPUT.in_file);
            ofstream ofs1("liberation1.xyz");
            ofstream ofs2("liberation2.xyz");
            ofstream ofs3("liberation3.xyz");
            ofstream ofs4("bend.xyz");
            ofstream ofs5("stretch.xyz");
            ofstream ofs6("astretch.xyz");
            ofstream ofs7("translation1.xyz");
            ofstream ofs8("translation2.xyz");
            ofstream ofs9("translation3.xyz");
            ofs1 << setprecision(8);
            ofs2 << setprecision(8);
            ofs3 << setprecision(8);
            ofs4 << setprecision(8);
            ofs5 << setprecision(8);
            ofs6 << setprecision(8);
            ofs7 << setprecision(8);
            ofs8 << setprecision(8);
            ofs9 << setprecision(8);
            
            string label;
            while(label != "ATOMIC_POSITIONS")
                Read_value(ifs,label); // find the position of atom positions
            
            Cellfile cel;
            Waterfile WF;
            
            cel.ReadInfile(ifs);
            WF.Read_water(cel);
            
            Waterana water(WF.waters[0],cel,INPUT.unitconv);
            Waterana tmpwater;
            double d = INPUT.displacement;
            
            //liberation1
            tmpwater = water.liberate(d);
            water.print_xyz(ofs1,"0d");
            tmpwater.print_xyz(ofs1,"d");
            tmpwater.liberate(d).print_xyz(ofs1,"2d");
            tmpwater.print_xyz(ofs1,"d");
            d *= -1;
            tmpwater = water.liberate(d);
            water.print_xyz(ofs1,"0d");
            tmpwater.print_xyz(ofs1,"-d");
            tmpwater.liberate(d).print_xyz(ofs1,"-2d");
            tmpwater.print_xyz(ofs1,"-d");
            d *= -1;
            
            //liberation2
            Vector3<double> axis;
            axis = water.OH[0] + water.OH[1];
            tmpwater = water.liberate(axis,d);
            water.print_xyz(ofs2,"0d");
            tmpwater.print_xyz(ofs2,"d");
            tmpwater.liberate(axis,d).print_xyz(ofs2,"2d");
            tmpwater.print_xyz(ofs2,"d");
            d *= -1;
            tmpwater = water.liberate(axis,d);
            water.print_xyz(ofs2,"0d");
            tmpwater.print_xyz(ofs2,"-d");
            tmpwater.liberate(axis,d).print_xyz(ofs2,"-2d");
            tmpwater.print_xyz(ofs2,"-d");
            d *= -1;
            
            //liberation3
            axis = axis.cross(water.vnorm);
            tmpwater = water.liberate(axis,d);
            water.print_xyz(ofs3,"0d");
            tmpwater.print_xyz(ofs3,"d");
            tmpwater.liberate(axis,d).print_xyz(ofs3,"2d");
            tmpwater.print_xyz(ofs3,"d");
            d *= -1;
            tmpwater = water.liberate(axis,d);
            water.print_xyz(ofs3,"0d");
            tmpwater.print_xyz(ofs3,"-d");
            tmpwater.liberate(axis,d).print_xyz(ofs3,"-2d");
            tmpwater.print_xyz(ofs3,"-d");
            d *= -1;
            
            //bending
            tmpwater.set(water.bend(d));
            water.print_xyz(ofs4,"0d");
            tmpwater.print_xyz(ofs4,"d");
            tmpwater.bend(d).print_xyz(ofs4,"2d");
            tmpwater.print_xyz(ofs4,"d");
            d *= -1;
            tmpwater.set(water.bend(d));
            water.print_xyz(ofs4,"0d");
            tmpwater.print_xyz(ofs4,"-d");
            tmpwater.bend(d).print_xyz(ofs4,"-2d");
            tmpwater.print_xyz(ofs4,"-d");
            d *= -1;
            
            //sstretching
            tmpwater.set(water.sstretch(d));
            water.print_xyz(ofs5,"0d");
            tmpwater.print_xyz(ofs5,"d");
            tmpwater.sstretch(d).print_xyz(ofs5,"2d");
            tmpwater.print_xyz(ofs5,"d");
            d *= -1;
            tmpwater.set(water.sstretch(d));
            water.print_xyz(ofs5,"0d");
            tmpwater.print_xyz(ofs5,"-d");
            tmpwater.sstretch(d).print_xyz(ofs5,"-2d");
            tmpwater.print_xyz(ofs5,"-d");
            d *= -1;
            
            //astretching
            tmpwater.set(water.astretch(d));
            water.print_xyz(ofs6,"0d");
            tmpwater.print_xyz(ofs6,"d");
            tmpwater.astretch(d).print_xyz(ofs6,"2d");
            tmpwater.print_xyz(ofs6,"d");
            d *= -1;
            tmpwater.set(water.astretch(d));
            water.print_xyz(ofs6,"0d");
            tmpwater.print_xyz(ofs6,"-d");
            tmpwater.astretch(d).print_xyz(ofs6,"-2d");
            tmpwater.print_xyz(ofs6,"-d");
            d *= -1;
            
            //translation1
            tmpwater = water.translate(water.vnorm*d);
            water.print_xyz(ofs1,"0d");
            tmpwater.print_xyz(ofs1,"d");
            tmpwater.translate(water.vnorm*d).print_xyz(ofs1,"2d");
            tmpwater.print_xyz(ofs1,"d");
            d *= -1;
            tmpwater = water.translate(water.vnorm*d);
            water.print_xyz(ofs1,"0d");
            tmpwater.print_xyz(ofs1,"-d");
            tmpwater.translate(water.vnorm*d).print_xyz(ofs1,"-2d");
            tmpwater.print_xyz(ofs1,"-d");
            d *= -1;
            
            //liberation;
            axis = water.OH[0] + water.OH[1];
            tmpwater = water.translate(axis*d);
            water.print_xyz(ofs2,"0d");
            tmpwater.print_xyz(ofs2,"d");
            tmpwater.translate(axis*d).print_xyz(ofs2,"2d");
            tmpwater.print_xyz(ofs2,"d");
            d *= -1;
            tmpwater = water.translate(axis*d);
            water.print_xyz(ofs2,"0d");
            tmpwater.print_xyz(ofs2,"-d");
            tmpwater.translate(axis*d).print_xyz(ofs2,"-2d");
            tmpwater.print_xyz(ofs2,"-d");
            d *= -1;
            
            //liberation3
            axis = axis.cross(water.vnorm);
            tmpwater = water.translate(axis*d);
            water.print_xyz(ofs3,"0d");
            tmpwater.print_xyz(ofs3,"d");
            tmpwater.translate(axis*d).print_xyz(ofs3,"2d");
            tmpwater.print_xyz(ofs3,"d");
            d *= -1;
            tmpwater = water.liberate(axis,d);
            water.print_xyz(ofs3,"0d");
            tmpwater.print_xyz(ofs3,"-d");
            tmpwater.translate(axis*d).print_xyz(ofs3,"-2d");
            tmpwater.print_xyz(ofs3,"-d");
            d *= -1;
        }
            break;
        case 2:
        {
            ifstream ifs(INPUT.in_file.c_str());\
            file_assert(ifs,INPUT.in_file);
            ofstream ofs("bending.xyz");
            ofs << setprecision(8);
            
            string label;
            while(label != "ATOMIC_POSITIONS")
                Read_value(ifs,label); // find the position of atom positions
            
            Cellfile cel;
            Waterfile WF;
            
            cel.ReadInfile(ifs);
            WF.Read_water(cel);
            Waterana water(WF.waters[0],cel,INPUT.unitconv);
    
            my_assert(INPUT.delta >=0, "Delta shoule >=0 !");
            //set mass of O and H to 16 and 2
            water.nbend(INPUT.delta,INPUT.displacement,INPUT.atom_mass[WF.O],INPUT.atom_mass[WF.H],true,ofs);
            
        }
            break;
        default:
            break;
    }
    
}

