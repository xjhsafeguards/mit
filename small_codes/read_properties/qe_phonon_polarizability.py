#!/usr/bin/env python3

import sys
import numpy as np

#read the polarizability as np array
def get_ndarray(line1,line2,line3):
    return np.array([line1.split()[3:6],line2.split()[3:6],line3.split()[3:6]],float)

ofs=open("alpha.eig.dat",'w')
ofs.write("#name pos1,2,3 distance alpha(1,2,3)\n")
for files in sys.argv:
    if(files.split('-',3)[-1]=='so-epsil-Ih-Cl-j'):
        files2=files.replace('so-epsil-Ih-Cl-j',files.split('-',3)[2]+'.xyz')
        print('reading files: ',files,' and ',files2)
        ifs1=open(files,'r')
        ifs2=open(files2,'r')
        line=ifs1.readline()
        #skip first line read the information to line2
        ifs2.readline()
        line2=ifs2.readline()
        ofs.write("#")
        ofs.write(line2)
        while(line):
            if 'Polarizability' in line :
                in_array=get_ndarray(ifs1.readline(),ifs1.readline(),ifs1.readline())
                #diagonalize the matrix
                w,v=np.linalg.eigh(in_array)
                #read the position
                pos=ifs2.readline().split()
                if(len(pos)>=5):
                    ofs.write('{} {} {} {} {} '.format(pos[0],pos[1],pos[2],pos[3],pos[4]))
                else:
                    ofs.write('{} {} {} {} {} '.format(pos[0],pos[1],pos[2],pos[3],0))
                ofs.write('{} {} {}\n'.format(w[0],w[1],w[2]))
            line=ifs1.readline()
        ifs1.close()
        ifs2.close()
ofs.close()







