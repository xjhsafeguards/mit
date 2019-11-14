/*************************************************************************
    > File Name: qe_internal_stress.cpp
    > Author: Jianhang Xu
    > Mail: jianhangxu@icloud.com 
    > Created Time: Thu Nov 14 12:34:32 2019
 ************************************************************************/

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cassert>
#include <ctime>

using namespace std;
template <class T>
static void READ_VALUE(ifstream &ifs, T &v)
{
	    ifs >> v;
	    ifs.ignore(150, '\n');
	    return;
}

int main()
{

	ifstream ifs_str("water.str");
	ifstream ifs_vel("water.vel");
	ifstream ifs_cel("water.cel");
	int natom = 192;
	int nss = 5; //105;

	int ss_index;
	double ss_time;
	double vx,vy,vz;

	double sx1, sx2, sx3;
	double sy1, sy2, sy3;
	double sz1, sz2, sz3;

	double cx1, cx2, cx3;
	double cy1, cy2, cy3;
	double cz1, cz2, cz3;

	double vol=0.0;
	double factor = 2.9421912*1.0e4*1822.8291744407; // a.u. to GPa

	ofstream ofs("out.ist");
	ofs << setprecision(6);

	for(int i=0; i<nss; ++i)
	{
		ifs_str >> ss_index >> ss_time;
		cout << ".str " << ss_index << " " << ss_time << endl;
		ifs_vel >> ss_index >> ss_time;
		cout << ".vel " << ss_index << " " << ss_time << endl;
		ifs_cel >> ss_index >> ss_time;
		cout << ".cel " << ss_index << " " << ss_time << endl;

		// unit is GPa
		ifs_str >> sx1 >> sx2 >> sx3;
		ifs_str >> sy1 >> sy2 >> sy3;
		ifs_str >> sz1 >> sz2 >> sz3;

		// unit is Bohr
		ifs_cel >> cx1 >> cx2 >> cx3;
		ifs_cel >> cy1 >> cy2 >> cy3;
		ifs_cel >> cz1 >> cz2 >> cz3;

		vol = cx1 * cy2 * cz3;
		cout << "volume is " << vol << " bohr^3" << endl;

		double kx1, kx2, kx3;
		double ky1, ky2, ky3;
		double kz1, kz2, kz3;
		kx1 = kx2 = kx3 = ky1 = ky2 = ky3 = kz1 = kz2 = kz3 = 0.0;
		for(int ia=0; ia<natom; ++ia)
		{
			double mass;
			if(ia<64) mass=15.9994;
			else mass=2.01410178;

			// I assume unit is in a.u. (bohr/a.u.)?
			// 1 a.u. = 2.9421912*10^13 Pa
			ifs_vel >> vx >> vy >> vz;
			kx1 += mass * vx * vx;
			kx2 += mass * vx * vy;
			kx3 += mass * vx * vz;
		
			ky1 += mass * vy * vx;
			ky2 += mass * vy * vy;
			ky3 += mass * vy * vz;
			
			kz1 += mass * vz * vx;
			kz2 += mass * vz * vy;
			kz3 += mass * vz * vz;
		}
		
		kx1 *= factor; kx2 *= factor; kx3 *= factor;
		ky1 *= factor; ky2 *= factor; ky3 *= factor;
		kz1 *= factor; kz2 *= factor; kz3 *= factor;

		kx1 /= vol; kx2 /= vol; kx3 /= vol;
		ky1 /= vol; ky2 /= vol; ky3 /= vol;
		kz1 /= vol; kz2 /= vol; kz3 /= vol;

		ofs << ss_index << " " << ss_time << endl;
		//ofs << kx1 << " " << kx2 << " " << kx3 << endl; 
		//ofs << ky1 << " " << ky2 << " " << ky3 << endl; 
		//ofs << kz1 << " " << kz2 << " " << kz3 << endl; 

		ofs << sx1-kx1 << " " << sx2-kx2 << " " << sx3-kx3 << endl; 
		ofs << sy1-ky1 << " " << sy2-ky2 << " " << sy3-ky3 << endl; 
		ofs << sz1-kz1 << " " << sz2-kz2 << " " << sz3-kz3 << endl; 
	}

	return 0;
}

