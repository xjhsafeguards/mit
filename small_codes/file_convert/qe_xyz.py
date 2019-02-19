#!/usr/local/bin/python3

nss=5
ssstep=0
natom=[64,128]
aname=['O','H']
unit_convert=0.529177249
inputfile='/Users/jianhangxu/Documents/2Cl/water_PBE/water64.pos'
outputfile='test.xyz'

def main():
    with open(inputfile,'r') as infile:
        with open(outputfile,'w') as outfile:
            for i in range(nss):
                pos2xyz(infile,outfile)
                for _ in range(ssstep)
                    next(infile)

def pos2xyz(ifs,ofs):
    totaln=0
    for nl in natom:
        totaln+=nl
    ofs.write('{}\n'.format(totaln))
    line = ifs.readline()
    ss = line.split()[0]
    time = line.split()[1]
    ofs.write('#' + ss + ' ' + time + '\n')
    for i in range(len(natom)):
        for j in range(natom[i]):
            ofs.write(aname[i])
            line = ifs.readline()
            for num in line.split():
                ofs.write('\t{}'.format(float(num)*unit_convert))
            ofs.write('\n')



if __name__ == '__main__':
    main()
