import numpy as np
from ase import io
from ase.units import *
import ase
import pickle

str2slice = lambda mystring: slice(*map(lambda x: int(x.strip()) if x.strip() else None, mystring.split(':')))

def read_lammps_dump(in_file="water.dump",natom=384,symbols="O128H256",nstep=10,nstart=1000,save_dump=True,rowcell=5,rowpos=9):
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

def read_ipi_xyz(infile="./",step=":",nbeads=8,unit=Bohr,colcell=2,save_dump=True,zfilldig=1):
    tmpq = []
    for i in range(0,nbeads):
        #print("Reading", i)
        tmpq.append(ase.io.read(infile+"data.pos_"+str(i).zfill(zfilldig)+".xyz",step))
    for beads in tmpq:
        for sn in beads:
            sn.set_pbc(111)
            #print(list(sn.info))
            #print(float(list(sn.info)[colcell]))
            cel = float(list(sn.info)[colcell])*unit
            sn.set_cell([cel,cel,cel])
            if(unit!=1):
                sn.set_positions(sn.get_positions()*unit)
    if(save_dump):
        pickle.dump(tmpq,open("data.pos.save","wb"))
    return tmpq

def read_dpmd_raw(infile="./",step=":",symbols="O64H128",unit=1,save_dump=True):
    c = []
    poss = np.loadtxt(infile+"coord.raw")
    boxs = np.loadtxt(infile+"box.raw")

    for pos,box in zip(poss[str2slice(step)],boxs[str2slice(step)]):
        #print(pos,box)
        c.append(ase.Atoms(symbols,cell=box.reshape(3,3)*unit,positions=pos.reshape(-1,3)*unit,pbc=(1,1,1)))
    if(save_dump):
        pickle.dump(c,open("raw.save","wb"))
    return c

from itertools import islice 

def read_cpmd(prefix="./water",step=":",symbols="O64H128",natoms=192,unit=Bohr,save_dump=True):
    c = []
    poss = open(prefix + ".pos","r").readlines()
    boxs = open(prefix + ".cel","r").readlines()
    
    ns = len(poss)//(natoms+1)
        
    pos_np_list = [np.array(line.split(),dtype=np.float) for line in poss]
    box_np_list = [np.array(line.split(),dtype=np.float) for line in boxs]
    
    itpos = iter(pos_np_list)
    itbox = iter(box_np_list)

    for pos,box in zip([list(islice(itpos,(natoms+1))) for _ in range(ns)][str2slice(step)],[list(islice(itbox,4)) for _ in range(ns)][str2slice(step)]):
        c.append(ase.Atoms(symbols,cell=np.array(box[1:])*unit,positions=np.array(pos[1:])*unit,pbc=(1,1,1)))
    if(save_dump):
         pickle.dump(c,open("cpmd.save","wb"))
    return c

#add by jianhang 06292020
import xml.etree.ElementTree as ET

def read_qboxr(infile="./",step=":",unit=Bohr,save_dump=True,if_pbc=True,speciesdict={"oxygen":"O","hydrogen":"H","chlorine":"Cl","sodium":"Na"}):
    tree = ET.parse(infile)
    root = tree.getroot()
    snap = []
    for it in root.findall("iteration")[str2slice(step)]:
        info=it.find("atomset")
        cell=info.find("unit_cell")
        in_cell = np.fromstring(cell.attrib["a"]+cell.attrib["b"]+cell.attrib["c"],sep=" ").reshape(3,3)*unit
        in_atoms = []
        for at in info.findall("atom"):
            in_atoms.append(ase.Atom(speciesdict.get(at.attrib["species"]),position=np.fromstring(at.find("position").text,sep=" ")*unit))
        snap.append(ase.Atoms(in_atoms,cell=in_cell,pbc=if_pbc))

    return snap

