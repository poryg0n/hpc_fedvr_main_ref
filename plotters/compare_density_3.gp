# Usage:
# gnuplot -c comp_density.gp datafile1 datafile2 datafile3 outputfile label1 label2

datafile1 = ARG1
datafile2 = ARG2
datafile3 = ARG3
outfile1  = ARG4
outfile2  = ARG5
label1    = ARG6 
label2    = ARG7 
label3    = ARG8 

#mode     = int(ARG3)

set format y "10^{%L}";

set terminal pngcairo size 1000,700
set output outfile1

set logscale y
set grid

set size ratio  0.3182 1,1
#set xrange [.5:1.13]
#set yrange [1.e-12:1.e0]
#set xtics .5, .1, 1.1
#set ytics 1.e-13, 1.e-2, 1.e-3

set xlabel 'x (a.u.)' font ', 15'
set ylabel '|{/Symbol y}(x,t)|^2' font ', 15'
set title  'Density Probability'

set grid

set key center bottom

plot datafile1.'/density.dat' using 1:2 with lines lw 2 lc 8 title label1, \
     datafile2.'/density.dat' using 1:2 with lines lw 2 lc 6 title label2, \
     datafile3.'/density.dat' using 1:2 with lines lw 2 lc 2 title label3


set terminal pngcairo size 1000,700
set output outfile2

set xlabel 'x (a.u.)' font ', 15'
#set ylabel 'D(x,t)' font ', 15'
set ylabel '|{/Symbol f}_{/Symbol w}(x,t)|^2' font ', 15'
set xlabel 'x (a.u.)'
set title  'Density Probability'

plot datafile1.'density.dat' using 1:3 with lines lw 2 lc 8 title label1, \
     datafile2.'density.dat' using 1:3 with lines lw 2 lc 6 title label2, \
     datafile3.'density.dat' using 1:3 with lines lw 2 lc 2 title label3 
