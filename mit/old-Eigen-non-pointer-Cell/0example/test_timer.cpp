#include <iomanip>
#include <chrono>
#include <cstring>
#include <functional>
#include <cmath>
#include <vector>

#include <Utility_time.h>

using namespace std;
using namespace chrono;


template<typename T>
void Print(const vector<T>& inv,ostream& os=cout){
    for( const auto& d: inv)
        //os << setw(10) << d;
        os << d << " ";
    os << endl;
}

double get_random(double min, double max) {
    /* Returns a random double between min and max */
    return (max-min) * ((double) rand() / (double) RAND_MAX) + min;
}

void MOD(double a,double b){
    cout << to_string(a) + " mod " + to_string(b) + " = " << remainder(a,b) << '\n';
}

void MOD2(double a,double b){
    cout << to_string(a) + " mod " + to_string(b) + " = " << fmod(a,b) << '\n';
}

double BC(double a, double b){
    while(a >= b/2){
        a -= b;
    }
    while(a < -b/2){
        a += b;
    }
    return a;
}

double BC2(double a, double b){
    while(a >= b){
        a -= b;
    }
    while(a < 0){
        a += b;
    }
    return a;
}

int main(int argc,char** argv){
    
    
    GTIMER.start("all");
    cout << "Modulus test\n";
    
    double a,b,r1,r2;
    
    for(int i=0; i< 1000000; ++i){
        //    GTIMER.start("random");
        a = get_random(-1000,1000);
        b = get_random(0,100);
        //    GTIMER.stop("random");
        //    GTIMER.start("cout");
        //cout << a << " " << b << '\r';
        //    GTIMER.stop("cout");
        //cout << "test1" << endl;
                GTIMER.start("reminder");
        //cout << "test2" << endl;
        r1 = remainder(a,b);
               GTIMER.stop("reminder");
        //cout << "test3" << endl;
                GTIMER.start("BC");
        //cout << "test4" << endl;
        r2 = BC(a,b);
                GTIMER.stop("BC");
        //cout << "test5" << endl;
        if((r1-r2)>0.00001)
            cout << a << " " << b << " " << r1 << " " << r2 << endl;
        //cout << "test6" << endl;
    }
    cout << endl;
    GTIMER.stop("all");
    GTIMER.summerize(cout);
    
}


//double cell_vol=0;
//int scount=0;
/*
 //auto Dp = new Distributionfunction(0,6,500);
 Distributionfunction* Dp[4];
 for(int i=0;i!=4;++i){
 //Dp[i] = new Distributionfunction(0,6,500);
 Dp[i] = new Distributionfunction(0,180,500); // OClO, ClOO, OOO, OO(cl-bonded)O
 }
 //for(int i=4;i!=6;++i){
 // Dp[i] = new Distributionfunction(0,180,500);
 //}
 
 for(int fc=0; fc!=8;++fc){
 
 ifstream ifs( file_folder + "/data.pos_"+ to_string(fc) + ".xyz");
 //molecule_manip* water = new water_manip();
 
 for(int i=0;i!=f_end;++i){
 
 std::shared_ptr<cell> cel = make_shared<cell_ipi>();
 cel->read(ifs);
 
 if(i>f_start and (i-f_start)%f_step == 0){
 
 vector<shared_ptr<position> > Os_in_Cl;
 for(int i=1;i<64;++i){
 bool is_in=false;
 if(cel->atoms()[0]->distance(*cel->atoms()[i])<OCl_cutoff){
 Os_in_Cl.push_back(cel->atoms()[i]);
 is_in=true;
 }
 vector<int> Oindex_in_Oi;
 for(int j=1;j<64;++j){
 if( i!=j and cel->atoms()[i]->distance(*cel->atoms()[j])<OO_cutoff)
 Oindex_in_Oi.push_back(j);
 }
 for(int j: Oindex_in_Oi){
 for(int k: Oindex_in_Oi)
 {
 if(j<k){
 double Angle=cel->atoms()[i]->angle(*cel->atoms()[j],*cel->atoms()[k]);
 Dp[2]->read(Angle); // OOO
 if(is_in)
 Dp[3]->read(Angle);
 }
 }
 }
 }
 for( const auto& atom1: Os_in_Cl){
 for( const auto& atom2: Os_in_Cl){
 if(atom1!=atom2 and atom1->distance(*atom2)<OO_cutoff)
 //data[0].push_back(cel->atoms()[0]->angle(*atom1,*atom2)); // OClO
 Dp[0]->read(cel->atoms()[0]->angle(*atom1,*atom2)); // OClO
 }
 for(int i=1;i<64;++i){
 if(atom1->distance(*cel->atoms()[i])<OO_cutoff)
 //data[1].push_back(atom1->angle(*cel->atoms()[0],*cel->atoms()[i])); // ClOO
 Dp[1]->read(atom1->angle(*cel->atoms()[0],*cel->atoms()[i])); // ClOO
 }
 }
 }
 cout << "read bead" << fc << " snapshot " << i << '\r' << flush;
 }
 }
 
 //cell_vol /= scount;
 
 vector<double> X=Dp[0]->get_x();
 vector<double> Y[4];
 for(int i=0;i!=4;++i){
 //Dp[i]->set_dimension(3);
 //Dp[i]->set_normalize(cell_vol);
 Y[i] = Dp[i]->get_y();
 }
 //for(int i=4;i!=6;++i){
 //Dp[i]->set_dimension(3);
 //Dp[i]->set_normalize(cell_vol);
 //Y[i] = Dp[i]->get_y();
 //}
 
 ofstream ofs("G_OClO_ClOO_OOO.txt");
 ofs << setprecision(10);
 
 ofs << "#" << setw(19) << "R" << setw(20) << "OClO" << setw(20) << "ClOO" << setw(20) << "OOO" << setw(20) << "OO(cl)O" << endl;
 
 for(int i=0;i<X.size();++i){
 ofs << setw(20) << X[i];
 for(int j=0;j!=4;++j){
 ofs << setw(20) << Y[j][i];
 }
 ofs << '\n';
 }
 */

