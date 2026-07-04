# Usage:
# gnuplot -c plot_pemd.gp datafile output1 output2 output3 logscale

datafile = ARG1
outfile1  = ARG2
outfile2  = ARG3
outfile3  = ARG4
logscale = int(ARG5)

set terminal pngcairo size 1000,700
set output outfile1

set xlabel 'k (a.u.)' font ', 15' enhanced
set ylabel 'P(k)' font ', 15' enhanced
set title 'Photoelectron Momentum Distribution (PEMD)' font ",16" enhanced

set format y "10^{%L}";

set grid
set size ratio 0.3182 1,1
set key left top
set xrange [-1.33:1.33]
set yrange [1.e-9 : 1.e-3]
set ytics  1.e-9, 1.e-2, 1.e-2
set xtics -1.0, .5, 1.0


if (logscale == 1) {
    set log y
#   plot datafile.'/pemd.dat' using 1:(log10($2)) with lines lw 2 lc 8 title 'P(k)'
    plot datafile.'/pemd.dat' using 1:($2) with lines lw 2 lc 8 title 'P(k)'
} else {
    plot datafile.'/pemd.dat' using 1:2 with lines lw 2 title 'P(k)'
}


set output outfile2
set ylabel 'b_{kw}(T)' font ', 15' enhanced
set yrange [1.e-6 : 1.e1]
set ytics  1.e-6, 1.e-3, 1.e3
#unset ytics
#unset format y

if (logscale == 1) {
    set log y
    plot datafile.'/pemd.dat' using 1:($4) with lines lw 2 lc 8 title '|b_{kw}(T)|^2', \
} else {
    plot datafile.'/pemd.dat' using 1:4 with lines lw 2 title '|b_{kw}(T)|^2'
}


set output outfile3
set ylabel 'Amplitudes' font ', 15' enhanced

    plot datafile.'/pemd.dat' using 1:($3*10**4) with lines lw 2 lc 8 title 'Re{a_k}', \
         datafile.'/pemd.dat' using 1:($4) with lines lw 2 lc 6 title 'Re{b_{kw}(T)}'
