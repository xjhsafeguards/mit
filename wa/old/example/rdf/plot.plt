#!/opt/local/bin/gnuplot

set terminal postscript color
set output "RDF_OO.eps"

set xrange [2:6]
set yrange [0:4]

set ylabel "g_{OO}(r)"
set xlabel "r[A]"
set grid

unset key
plot "RDF_OO.txt" using 1:2 with l lw 3,

set output "RDF_OH.eps"

set xrange [0.5:4]
set yrange [0:2]

set ylabel "g_{OH}(r)"
set xlabel "r[A]"
set grid

unset key
plot "RDF_OH.txt" using 1:2 with l lw 3,
unset output
