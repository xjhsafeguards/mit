
awk 'BEGIN {for(i=0;i<128;i++) printf "0 "; for(i=0;i<256;i++) printf "1 " }' > type.raw

INFILE=scan.out
force=25.7105
energy=13.6056925
length=0.52917721092
box_bohr=29.58541
box_angstrom1=$box_bohr*$length
box_angstrom2=$box_angstrom1
box_angstrom3=$box_angstrom1

#awk 'BEGIN {ORS=" "} {printf "%.10f %.10f %.10f ",$3,$4,$5}' test-cla.xyz > position.raw
awk 'BEGIN {ORS=" "} {printf "%.10f %.10f %.10f ",$2*'$length',$3*'$length',$4*'$length'}' cell.xyz > coord.raw


#####
grep atom $INFILE |grep force | awk 'BEGIN {ORS=" "} {printf "%.10f %.10f %.10f ",$7*'$force',$8*'$force',$9*'$force'}' > force.raw
grep ! $INFILE| tail -n1|awk '{printf "%.10f ",$5*'$energy'}' > energy.raw
echo | awk '{printf "%.10f 0 0 0 %.10f 0 0 0 %.10f ",'$box_angstrom1','$box_angstrom2','$box_angstrom3'}' > box.raw
