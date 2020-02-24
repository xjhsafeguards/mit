import numpy as np
from ase import io
from ase.units import *
import ase
import pickle

def readi_lammps_dump(in_file="water.dump",natom=384,symbols="O128H256",nstep=10,nstart=1000,save_dump=True,rowcell=5,rowpos=9):
    c = []
    dumpfile = open(in_file,"r")
    content  = dumpfile.readlines()
    lines = len(content)
    ns = lines // (natom+rowpos)

    for i in range(ns):
        if(i>nstart and i%nstep == 0):
            snapstart = i*(natom+rowpos)
            cell = np.array(content[snapstart+rowcell].split()[:2],dtype=np.float)
            cellsize = cell[1] - cell[0]
            pos = np.zeros((natom,5))
            for j in range(natom):
                pos[j]=np.array(content[snapstart+rowpos+j].split(),dtype=np.float)
            c.append(ase.Atoms(symbols,cell=[cellsize,cellsize,cellsize],positions=pos[pos[:,0].argsort(),-3:],pbc=(1,1,1)))
    if(save_dump):
        pickle.dump(c,open(in_file+".save","wb"))
    return c

def read_ipi_xyz(infile="./",step=":",nbeads=8,unit=Bohr,colcell=3,save_dump=True):
    tmpq = []
    for i in range(0,nbeads):
        #print("Reading", i)
        tmpq.append(ase.io.read(infile+"data.pos_"+str(i)+".xyz",step))
    for beads in tmpq:
        for sn in beads:
            sn.set_pbc(111)
            cel = list(sn.info)[colcell]*unit
            sn.set_cell([cel,cel,cel])
            if(unit!=1):
                sn.set_positions(sn.get_positions()*unit)
    if(save_dump):
        pickle.dump(tmpq,open("data.pos.save","wb"))
    return tmpq
