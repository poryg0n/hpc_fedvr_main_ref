# Usage:
# gnuplot -c plot_density.gp datafile1 datafile2 outfile1 outfile2 label1 label2


datafile1 = ARG1
datafile2 = ARG2
outfile1  = ARG3
outfile2  = ARG4
label1  = ARG5
label2  = ARG6
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
set title  'Density Probability'

set grid

set key center bottom

#plot for [i=1:ARGC] \
#     ARG1/density.dat using 1:2 with lines lw 2 lc 8 title "|{/Symbol y}(x,T)|^2", \
#     ARG2/density.dat using 1:2 with lines lw 2 lc 6 title "|{/Symbol y}(x,T)|^2", \
#    datafile using 1:3 with lines lw 2 lc 2 title "|{/Symbol f}_{/Symbol w}(x,T)|^2"


plot datafile1 using 1:2 with lines lw 2 lc 8 title label1, \
     datafile2 using 1:2 with lines lw 2 dt 2 lc 6 title label2, \


set output outfile2
plot datafile1 using 1:3 with lines lw 2 lc 8 title "|{/Symbol f}_{/Symbol w}(x,T)|^2", \
     datafile2 using 1:3 with lines lw 2 dt 2 lc 6 title "|{/Symbol f}_{/Symbol w}(x,T)|^2", \
