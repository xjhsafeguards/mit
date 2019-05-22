#include "Base_cell.h"

using namespace std;

int main(){
    Box b1,b2;
    //cout << "input your cell" << endl;
    //b1.read(cin);
    cout << "your cell b2:" << endl;
    //b1.write(cout) << endl;
    b2.read(3,4,3);
    b2 << 1,2,3,4,5,6,7,8,10;
    cout << b2.inverse() << endl;
    b2.write(cout) << endl;
    b2 *=2;
    b2.write(cout) << endl;
    cout << b2.volume() << endl;
    Vector v1,v2,v3;
    v1 << 0,0,0;
    v2 << 0.5,0,0;
    v3 << 0,0.5,0;
    cout << "v1,v2,v3:\n";
    cout << v1 << '\n' << v2 << '\n' << v3 << endl;
    cout << b2.shortest_fvector(v1,v2) << endl;
    cout << b2.cdistance(v1,v2) << " " << b2.fdistance(v1,v2) << endl;
    cout << b2.cangle(v1,v2,v3) << " " << b2.fangle(v1,v2,v3) << endl;
    cout << b2.cangle(v2,v1,v3) << " " << b2.fangle(v2,v1,v3) << endl;
    cout << b2.to_fposition(b2.to_cposition(b2.to_fposition(b2.to_cposition(v2))))<< endl;
    
    {
        cout << "test for positions" << endl;
        Positions p1;
        p1.init(2);
        p1 << 1,2,3,4,5,6;
        cout << p1 << endl;
        p1.read(1,2,3,4);
        p1.write(cout)<< endl;
        //cout << "put your first line of pi: ";
        //p1.read(0,cin);
        //p1.write(cout) << endl;
        //cout << p1(1) << endl;
        
        const Vector& r1= p1.row(0),& r2=p1.row(1);
        //r1 << 0,1,2;
        cout << r1 << endl << r2 << endl;
        cout << p1 << endl;
    }
    
    {
        cout << "test for position" << endl;
        Position p1(0,0,0),p2;
        p2 << 4,5,6;
        p2.read(1,1,1);
        cout << p1 << p2 << endl;
        cout << b2 << endl;
        cout << b2.cdistance(p1,p2) << " " << b2.fdistance(p1,p2) << endl;
        Position_pv pv1;
        pv1.add(0,0,0).add(cin);
        cout << *pv1[1] << endl;
    }

}
