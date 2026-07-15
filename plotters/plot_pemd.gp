# Usage:
# gnuplot -c plot_pemd.gp datafile output1 output2 output3 logscale

datafile = ARG1
outfile1  = ARG2
outfile2  = ARG3
outfile3  = ARG4

set terminal pngcairo size 1000,700
set output outfile1

set xlabel 'k (a.u.)' font ', 15' enhanced
set ylabel 'P(k)' font ', 15' enhanced
set title 'Photoelectron Momentum Distribution (PEMD)' font ",16" enhanced

set format y "10^{%L}";
set log y

set grid
set size ratio 0.3182 1,1
set key left top
set xrange [-1.33:1.33]
set yrange [1.e-9 : 1.e-3]
set ytics  1.e-9, 1.e-2, 1.e-2
set xtics -1.0, .5, 1.0


    plot datafile.'/pemd.dat' using 1:2 with lines lw 2 lc 8 title 'P(k)'


set output outfile2
set ylabel 'b_{kw}(T)' font ', 15' enhanced
set yrange [1.e-6 : 1.e1]
set ytics  1.e-6, 1.e-3, 1.e3
#unset ytics
#unset format y

set output outfile2
set ylabel 'Amplitudes' font ', 15' enhanced
set key center

    plot datafile.'/pemd.dat' using 1:($2*10**3) with lines lw 2 lc 8 title '|a_k|^2', \
         datafile.'/pemd.dat' using 1:($3) with lines lw 2 lc 6 title '|b_{kw}(T)|^2'
