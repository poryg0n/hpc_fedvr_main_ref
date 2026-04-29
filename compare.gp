# Usage:
# gnuplot -c plot_density.gp datafile outputfile

datafile1 = ARG1
datafile2 = ARG2
outfile1  = ARG3
outfile2  = ARG4
label1    = ARG5 
label2    = ARG6 

#mode     = int(ARG3)

set format y "10^{%L}";

set terminal pngcairo size 1000,700
set output outfile1

set logscale y
set grid

set size ratio  0.3182 1,1
#set xrange [.5:1.13]
#set yrange [1.e-14:1.e-2]
#set xtics .5, .1, 1.1
#set ytics 1.e-13, 1.e-2, 1.e-3

set xlabel 'x (a.u.)' font ', 15'
set ylabel 'D(x,t)' font ', 15'
set xlabel 'x (a.u.)'
set ylabel 'D(x,T)'
set title  'Density Probability'

set grid

set key center bottom

plot datafile1 using 1:2 with lines lw 2 lc 8 title label1, \
     datafile2 using 1:2 with lines lw 1 lc 2 title label2


set terminal pngcairo size 1000,700
set output outfile2

set xlabel 'x (a.u.)' font ', 15'
set ylabel 'D(x,t)' font ', 15'
set xlabel 'x (a.u.)'
set ylabel 'D(x,T)'
set title  'Density Probability'

plot datafile1 using 1:3 with lines lw 2 lc 8 title label1, \
     datafile2 using 1:3 with lines lw 1 lc 2 title label2 
