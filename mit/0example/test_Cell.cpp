#include "Cell.h"

using namespace std;

int main(){
    Box b1;
    //cout << "please input cell" << endl;
    b1.load(3,1,1,2,6,3,1,4,9);
    cout << "box input: \n" << b1 << endl;
    cout << "volume of the box: " << b1.volume() << endl;
    Eigen::Matrix<double,1,3,Eigen::RowMajor> p1,p2,p3;
    Positions ps1;
    Cell cel1;
    ps1.resize(3,3);
    p1 <<    0.991875736,0.849504039,0.994397779;
    p2 <<    0.998704088,0.933892637,0.884919008;
    p3 <<    0.086486646,0.161505974,0.852384876;
    ps1.row(0) << p1;
    ps1.row(1) << p2;
    ps1.row(2) << p3;
    cout << b1.angle(p1,p2,p3) << " " << b1.distance(p1,p2) << endl;
    for(int i=0;i<3;++i){
        cout << b1.distance(ps1.pos(i),p1) << endl;
    }
    cel1.cellp = b1;
    cel1.atoms = ps1;
    ps1.row(0) << 0,0,0;
    for(int i=0;i<3;++i){
        cout << cel1.distance(cel1.atom(i),p1) << endl;
    }
    for(int i=0;i<3;++i){
        cout << ps1.vpos(0)[i] << " ";
    }
    cout << endl << ps1.pos(0) << endl;
    for(int i=0;i<3;++i){
        cout << ps1.vpos(1)[i] << " ";
    }
    cout << endl << ps1.pos(1) << endl;
    
    cel1.atoms.cart(cel1.cellp);
    cout << cel1.atoms << endl;
    cel1.atoms.frac(cel1.cellp);
    cout << cel1.atoms << endl;
}
