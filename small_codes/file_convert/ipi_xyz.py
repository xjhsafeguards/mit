#!/usr/local/bin/python3

#for qe
nss=1937
ssstep=10
natom=[64,128]
aname=['O','H']

#for ipi

#global
unit_convert=0.529177249
inputfile='data.pos_0.xyz'
outputfile='test.xyz'

def main():
    with open(inputfile,'r') as infile:
        with open(outputfile,'w') as outfile:
            for i in range(nss):
                ipi2xyz(infile,outfile)


def ipi2xyz(ifs,ofs):
    line = ifs.readline()
    ipi_natom = int(line.split()[0])
    ofs.write('{}\n'.format(ipi_natom))
    line = ifs.readline()
    ofs.write(line)
    for i in range(ipi_natom):
        line = ifs.readline()
        ofs.write(line.split()[0])
        for num in line.split()[1:]:
            ofs.write('\t{:12.8f}'.format(float(num)*unit_convert))
        ofs.write('\n')
    for _ in range(ssstep):
        for _ in range(ipi_natom+2):
            next(ifs)


if __name__ == '__main__':
    main()
