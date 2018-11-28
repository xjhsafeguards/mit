#!/opt/local/bin/gnuplot

set terminal postscript color
set output "IR.eps"

#set xrange [2000:2800]

set ylabel "a(w)n(w)[10^3 cm^{-1}]"
set xlabel "wavenumber[cm^{-1}]"
set grid

unset key
plot "FT.txt" using ($1):2 with l lw 3,

unset output
