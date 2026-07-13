# Usage:
# gnuplot -c plot_density.gp datafile outputfile


datafile = ARG1
outfile  = ARG2
#mode     = int(ARG3)

#set format y "10^{%L}";

set terminal pngcairo size 1000,700
set output outfile

#set logscale y
set grid

set size ratio  0.3182 1,1
set xrange [-10.:10.]
#set yrange [1.e-14:1.e-2]
#set xtics .5, .1, 1.1
#set ytics 1.e-13, 1.e-2, 1.e-3

set xlabel 'x (a.u.)' font ', 15'
set ylabel 'wf' font ', 15'
set title  'Wavefunction'

set grid

set key center bottom

plot datafile."init_psi.dat" using 1:2 with linesp dt 2 lw 2 lc 8 title "{/Symbol y}(x,T)", \
     datafile."init_phi.dat" using 1:3 with linesp dt 2 lw 2 lc 2 title "{/Symbol f}_{/Symbol w}(x,T)"
