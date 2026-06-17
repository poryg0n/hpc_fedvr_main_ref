# Usage:
# gnuplot -c plot_pemd.gp datafile output1 output2 logscale

datafile = ARG1
outfile1  = ARG2
#outfile2  = ARG3
#logscale = int(ARG4)

set terminal pngcairo size 1000,700
set output outfile1

set xlabel 't (a.u.)' font ', 15' enhanced
set ylabel '{/Symbol \362} D(x,t)dx' font ', 15' enhanced
set ylabel '{/Symbol \326}{/Symbol S} D(x,t)' font ', 15' enhanced
set ylabel 'norm' font ', 15' enhanced
#set ylabel '$ D(x,t)dx$' font ', 15' enhanced
set title 'Evolution of the norm' font ",16" enhanced

#set format y "10^{%L}";

set grid
set size ratio 0.3182 1,1
set key left top
#set xrange [-1.33:1.33]
#set yrange [1.e-14 : 1.e-3]
#set ytics  1.e-9, 1.e-2, 1.e-2
#set xtics -1.0, .5, 1.0

#set log y
plot datafile."/norm.dat" using 1:($2) with lines lw 2 lc 6 title "|{/Symbol y}(x,t)|", \
     datafile."/norm.dat" using 1:($3) with lines lw 2 lc 8 title "|{/Symbol f}_{/Symbol w}(x,t)|", \
