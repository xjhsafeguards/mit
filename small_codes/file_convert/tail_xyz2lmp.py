import numpy as np
pos = np.loadtxt("pos_angs.xyz")
from Cl_head import *
pos_angs = pos*Bohr
col1 = np.arange(190) + 1
col2 = np.concatenate((np.full(1,1),np.full(63,2),np.full(126,3)))
with open('tail.lmp','w') as ofs:
    for i in range(190):
        ofs.write(str(col1[i]))
        ofs.write("     ")
        ofs.write(str(col2[i]))
        ofs.write("     ")
        ofs.write(str(pos_angs[i])[1:-1])
        ofs.write("\n")
