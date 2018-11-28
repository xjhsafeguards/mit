#!/gnuplot

set terminal postscript color 

set output 'dens_z.eps'
#set style fill solid 0.4
set xlabel 'z-direction(A)'
set ylabel 'density(g/cm^3)'
#set yrange [0:1]
set title 'Density Distribution'
#unset key

plot 'density_Gs.txt' u 1:2 w l

unset output
