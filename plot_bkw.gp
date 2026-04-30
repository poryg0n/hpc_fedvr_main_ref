# Usage:
# gnuplot -c plot_bkw.gp datafile outputfile

datafile = ARG1
outfile  = ARG2
#mode     = int(ARG3)

set format y "10^{%L}";

set terminal pngcairo size 1000,700
set size ratio  0.3182 1,1

set output outfile

set datafile commentschars "#"


set title "Momentum spectrum"
set xlabel "k (a.u.)"
set ylabel "|b_{k{/Symbol w}}|^2"
set xrange [-1.6:1.6]

set grid

set logscale y
set format y "10^{%L}"

plot \
    datafile using 1:(($2)**2 + ($3)**2) w l lw 2 title "|b_k|^2"
