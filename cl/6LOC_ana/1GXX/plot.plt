#!/opt/local/bin/gnuplot

set terminal postscript enhanced color


set output "Cl-O.eps"

set grid
set grid mxtics mytics
set mxtics 5 
set mytics 5

set style line 1 lw 3

set xrange [2:6]

QUANTUM="../../../nvt_cl_63H2O_pi_300K/data_720/GXX"
CLASSICAL="."
EXPR="/Users/jianhangxu/Library/Mobile Documents/com~apple~CloudDocs/study/work/WCl/1structure/2G_XX/expr"

plot QUANTUM."/G_OH_OO_OCl_HCl.txt" using 1:4 with l lw 3 t "quantum",\
CLASSICAL."/G_OH_OO_OCl_HCl.txt" using 1:4 with l lw 3 t "classical",\
EXPR."/G_Cl_O.txt" using 1:2 with l lw 3 t "expr"


set output "Cl-H.eps"

set xrange [1:6]

plot QUANTUM."/G_OH_OO_OCl_HCl.txt" using 1:5 with l lw 3 t "quantum",\
CLASSICAL."/G_OH_OO_OCl_HCl.txt" using 1:5 with l lw 3 t "classical",\
EXPR."/G_Cl_H.txt" using 1:2 with l lw 3 t "expr"

set output "O-O.eps"

set xrange [2:6]

plot QUANTUM."/G_OH_OO_OCl_HCl.txt" using 1:3 with l lw 3 t "quantum",\
CLASSICAL."/G_OH_OO_OCl_HCl.txt" using 1:3 with l lw 3 t "classical",


set output "O-H.eps"

set xrange [0:6]
set yrange [0:3]

plot QUANTUM."/G_OH_OO_OCl_HCl.txt" using 1:2 with l lw 3 t "quantum",\
CLASSICAL."/G_OH_OO_OCl_HCl.txt" using 1:2 with l lw 3 t "classical",

set output "O-Cl-O.eps"

set xrange [0:180]
set yrange [0:0.02]

plot QUANTUM."/G_OH_OO_OCl_HCl.txt" using ($1*30):6 with l lw 3 t "quantum",\
CLASSICAL."/G_OH_OO_OCl_HCl.txt" using ($1*30):6 with l lw 3 t "classical",\

set output "H-Cl-H.eps"

set xrange [0:180]
set yrange [0:0.02]

plot QUANTUM."/G_OH_OO_OCl_HCl.txt" using ($1*30):7 with l lw 3 t "quantum",\
CLASSICAL."/G_OH_OO_OCl_HCl.txt" using ($1*30):7 with l lw 3 t "classical",\
EXPR."/P_H-Cl-H.txt" using ($1):2 with l lw 3 t "exp",

set output "H-Cl-H_comp.eps"

set xrange [0:180]
set yrange [0:0.02]

plot QUANTUM."/G_OH_OO_OCl_HCl.txt" using ($1*30):7 with l lw 3 t "SCAN quantum",\
CLASSICAL."/G_OH_OO_OCl_HCl.txt" using ($1*30):7 with l lw 3 t "SCAN classical",\
"/Users/jianhangxu/Library/Mobile Documents/com~apple~CloudDocs/study/work/WCl/PBE/3H-Cl-H/H_cutoff/quantum_2.0_4.8/adf_Hcutoff.txt" using 1:7 with l lw 3 t "PBE quantum",\
"/Users/jianhangxu/Library/Mobile Documents/com~apple~CloudDocs/study/work/WCl/PBE/3H-Cl-H/H_cutoff/classical/adf_Hcutoff.txt" using 1:7 with l lw 3 t "PBE classical",\
EXPR."/P_H-Cl-H.txt" using ($1):2 with l lw 3 t "exp",

