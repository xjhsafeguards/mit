#!/usr/local/bin/python3

nss=100000
ssstep=0
natom=[160,320]
aname=['O','H']
unit_convert=0.529177249
inputfile='/Users/jianhangxu/Documents/My_cpmd/H2O-vacum/water-vacuum-scan-nvt-300K-160/H2O.pos'
cellfile='/Users/jianhangxu/Documents/My_cpmd/H2O-vacum/water-vacuum-scan-nvt-300K-160/H2O.cel'
outputfile='test.xyz'

def main():
    with open(cellfile,'r') as cfile:
        with open(inputfile,'r') as infile:
            with open(outputfile,'w') as outfile:
                totaln=0
                for nl in natom:
                    totaln+=nl
                for i in range(nss):
                    pos2xyz(infile,cfile,outfile)
                    for _ in range(ssstep):
                        for _ in range(totaln+1):
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

def pos2xyz(ifs,ifs1,ofs):
    totaln=0
    for nl in natom:
        totaln+=nl
    ofs.write('{}\n'.format(totaln))
    line = ifs.readline()
    line1 = ifs1.readline()
    ss = line.split()[0]
    ss1 = line1.split()[0]
    if(ss!=ss1):
        os.exit("pos cel dont match in ss" + str(ss) + " !")
    time = line.split()[1]
    ofs.write('#' + ss + ' ' + time + '\n')
    cell = []
    cell.append(ifs1.readline().split()[0])
    cell.append(ifs1.readline().split()[1])
    cell.append(ifs1.readline().split()[2])
#print(cell)
    for i in range(len(natom)):
        for j in range(natom[i]):
            ofs.write(aname[i])
            line = ifs.readline()
            for num in range(len(line.split())):
                ofs.write('\t{}'.format(BC(float(line.split()[num]),float(cell[num]))*unit_convert))
            ofs.write('\n')

def BC(num,bc):
    result=num
    while(result>=bc):
        result=result-bc
    #print(result)
    while(result<0.0):
        result=result+bc
    #print(result)
    return result


if __name__ == '__main__':
    main()
