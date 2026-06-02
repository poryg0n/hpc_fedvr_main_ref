# gnuplot -c "plot_pk0.gp" datafile outfile
# Usage:

datafile1 = ARG1
datafile2 = ARG2
outfile1  = ARG3
outfile2  = ARG4
label1    = ARG5
label2    = ARG6


set terminal pngcairo size 1000,700

#mode     = int(ARG3)

set terminal pngcairo size 1000,700
set output outfile1

set title "p_{k0} convergence check"
set xlabel "k (a.u.)" font ", 15"
set ylabel "p_{k0}" font ", 15"

set size ratio  0.3182 1,1
#set xrange [.5:1.13]
#set yrange [1.e-14:1.e-2]
#set xtics .5, .1, 1.1
#set ytics 1.e-13, 1.e-2, 1.e-3


set grid
set key left top

plot \
datafile1 using 1:4 with lines lw 2 lc 8 title label1, \
datafile2 using 1:4 with lines lw 2 lc 6 title label2, \


set output outfile2

plot \
datafile1 using 1:2 with lines lw 2 lc 8 title label1, \
datafile2 using 1:3 with lines lw 2 lc 6 title label2, \
