#include "ghelp.h"

map<string,vector<string> > Help::dic = {
    {
        "hbs",
        {
            "Analyzing H_bond ",
            "list water and ion H-bonds distribution ",
            "list only water H-bonds ",
            "list only water H-bonds in z range INPUT.parameter[0,1] (angstrom)"
        }
    },
    {
        "wair",
        {
            "Analyzing IR spectrum of water ",
            "calculate general IR spectrum",
            "calculate spacial IR spectrum in INPUT.parameter[0,1] (Angstrom)",
            "calculate spacial IR spectrum in INPUT.parameter[0,1] (Angstrom) with given number INPUT.index hbond connection to first shell of the center",
            "calculate spacial IR spectrum in INPUT.parameter[0,1] (Angstrom) with given number INPUT.index hbond connection to center ",
            "calculate IR spectrum with given number INPUT.index hbond connection",
        }
    },
    {
        "dftcell",
        {
            "Analyzing Cell files "
            "Convert files from in_file to out_file_type ",
        }
    }
};

void Help::Routine(string calculation){
    if(Help::dic.count(calculation)){
        Help::title(calculation);
        for(int i=0;i<Help::size(calculation);i++)
            Help::content(calculation,i);
    }else
        cout << "No help section for calculation " << calculation << endl;
}

void Help::help(string calculation,int method){
    if(Help::dic.count(calculation) and Help::size(calculation) > method){
        Help::title(calculation);
        Help::content(calculation,method);
    }else
        cout << "Wrong type for calculation " << calculation << endl;
}

void Help::title(string calculation){
    cout << Help::dic[calculation][0] << "with type :" << endl;
}

void Help::content(string calculation, int method){
    cout << " ** " << method << " ** "<< Help::dic[calculation][method+1]<< endl;
}

int  Help::size(string calculation){
    return Help::dic[calculation].size() - 1;
}
/*
if(calculation=="density")
{
    cout << "Compute Density of the System with method: " << endl;
    cout << "** 0 ** Calculate both GBS and ILI" << endl;
    cout << "** 1 ** Calculate only the z axis density distribution" << endl;
    cout << "** 2 ** Calculate density distribution with GBS" << endl;
    cout << "** 3 ** Caluclate density distribution with ILI" << endl;
    
    }
    else if(calculation=="reorganize")
    {
        cout << "Reorganizing the System with type: " << endl;
        cout << "** 0 ** Generate .xyz file normally!" << endl;
        cout << "** 1 ** Generate .xyz file with non-water highlighted!" << endl;
        cout << "** 2 ** Generate .xyz file with non-water only!" << endl;
        cout << "** 3 ** Generate the compatiable .pos and .wfc files" << endl;
        cout << "** 4 ** Generate the average value of NPT cell" << endl;
    }
    else if(calculation=="hbs")
    {
        cout << "Analyzing H_bond using type: " << endl;
        cout << "** 0 ** list ion H-bonds distribution" << endl;
        cout << "** 1 ** list only water H-bonds" << endl;
    }
    else if(calculation=="ir")
    {
        cout << "Using Wannier Center Calculating dipole using type: " << endl;
        cout << "** 0 ** Start from Reading .wfc .pos (.vel) files" << endl;
        cout << "** 1 ** Generate dipole.txt and vdipole.txt only" << endl;
        cout << "** 2 ** Start from Reading vdipole.txt file" << endl;
        cout << "** 3 ** Calculate Intramolecular contribution from vdipole.txt file" << endl;
        cout << "** 4 ** Calculate electronic and ionic contributions" << endl;
        cout << "** 41 ** Calculate lone-pair and bond-pair electronic contributions" << endl;
        cout << "** 5 ** Generate electronic and ionic contributions to edipole.txt and idipole.txt only" << endl;
        cout << "** 51 ** Generate lone-pair and bond-pair electronic contributions to eldipole.txt and ebdipole.txt only" << endl;
        cout << "** 6 ** Calculate electronic and ionic contributions from edipole.txt and idipole.txt" << endl;
        cout << "** 7 ** Calculate elctronic ionic coupling term from ion,electron_vdipole.txt" << endl;
        cout << "** 71 ** Calculate lone-pair and bond-pair electronic coupling term from electron_lone/bond_vdipole.txt" << endl;
    }
    else if(calculation=="print_water")
    {
        cout << "Printing water of the System using type: "<< endl;
    }
    else if(calculation=="constrain_water")
    {
        cout << "Adding constrain on water to the System" << endl;
    }
    else if(calculation=="animate_water")
    {
        cout << "Animate water using type: "<< endl;
        cout << "Making water vibration eigen mode animation" << endl;
        cout << "Bending water with displacement and delta: " << endl;
    }
    else if(calculation=="rdf")
    {
        cout << "Calculating RDF using type: " << endl;
        cout << "** 0 ** Reading RDF for OO and OH in 1D" << endl;
        cout << "** 1 ** Reading RDF for OO and OH in 3D" << endl;
        cout << "** 2 ** Reading RDF for OMLWF in 1D" << endl;
        cout << "** 3 ** Reading RDF for OMLWF in 3D" << endl;
    }
    else if(calculation=="count")
    {
        cout << "Count how many snapshots were included in .wfc .pos .cel files." << endl;
    }
    else if(calculation=="hba")
    {
        cout << "hba is not fixed yet !" << endl;
    }
    else if(calculation=="wfa")
    {
        cout << " Analyze  water file using type: " << endl;
    }
    else if(calculation=="wfat")
    {
        cout << " Analyze time dependent water file using type: " << endl;
        cout << "** 0 ** Calculation average 2nd order Rotational Correlation Time of the O-H covalent bond! " << endl;
        cout << "** 1 ** Calculation 2nd order Rotational Correlation Time of the O-H covalent bond of certain O and H "<< endl;
    }
    else
    {
        cout << "No help section for calculation " << calculation << endl;
        }
    */
