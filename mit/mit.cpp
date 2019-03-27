#include "Cell.h"

#include <iostream>
#include <string>

using namespace std;

int main(int argc, char** argv){
    Cell cel;
    cel.cell_parameters = Eigen::Vector3d(10,10,10).asDiagonal();
    cel.positions.resize(3,3);
    cel.positions << 0,0,0, 3,0,0 ,3,4,0;
    cel.is_frac=false;
    
    cout << "cart: " << endl;
    cout << cel.positions << endl;
    
    cel.to_frac();
    cout << "frac: " << endl;
    cout << cel.positions << endl;
    
    cel.to_cart();
    cout << "cart: " << endl;
    cout << cel.positions << endl;
    
    cel.to_frac();
    cout << "frac: " << endl;
    cout << cel.positions << endl;
    
    cout << "angle 0,1,2:";
    cout << cel.positions.row(0) << endl;
    cout << cel.cal_angle(cel.positions.row(0),cel.positions.row(1),cel.positions.row(2),cel.cell_parameters) << endl;
    
    cout << "distance:";
    cout << cel.cal_distance(cel.positions.row(0),cel.positions.row(2),cel.cell_parameters) << endl;
    return 1;
}
