#!/usr/bin/env python3

import numpy as np
from io import StringIO

filename = "data.pos_1.xyz"
infile=open(filename,'r')

while(infile.readline()!=''):
    infile.readline()
#m = np.loadtxt(infile,dtype=[('type','S10'),('x','f'),('y','f'),('z','f')],max_rows=190)
    m = np.loadtxt(infile,usecols=(1,2,3),max_rows=190)
#print(m)
    #for line in open(filename,'r'):
    #m = np.loadtxt(StringIO(line),dtype=None)
'''with open(filename,'r') as infile:
    for line in infile:
        pass
        '''
#data = np.genfromtxt(filename)
