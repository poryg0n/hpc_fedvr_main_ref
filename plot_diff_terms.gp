# Usage:
# gnuplot -c plot_pemd.gp datafile output1 output2 logscale

datafile = ARG1
outfile1  = ARG2
outfile2  = ARG3
logscale = int(ARG4)

set terminal pngcairo size 1000,700
set output outfile1

set xlabel 'k (a.u.)' font ', 15' enhanced
set ylabel 'C_{k{/Symbol w}} p_{k0}a_0' font ', 15' enhanced
set title 'Non integral term in b_{k{/Symbol w}}' font ",16" enhanced

set format y "10^{%L}";

set grid
set size ratio 0.3182 1,1
set key left top
#set xrange [-1.33:1.33]
#set yrange [1.e-14 : 1.e-3]
#set ytics  1.e-9, 1.e-2, 1.e-2
#set xtics -1.0, .5, 1.0


#if (logscale == 1) {
#    set log y
## plot the non integral constant term
#    plot datafile using 1:($2) with lines lw 2 lc 7 title 'C_{k{/Symbol w}}p_{k0}a_0'
#} else {
#    plot datafile using 1:($2) with lines lw 2 lc 7 title 'C_{/Symbol w}_{k0}p_{k0}a_0'
#}

#
set log y
plot datafile."/vec_01k.dat" using 1:(abs($2)) with lines lw 2 lc 7 title "Re(Int p_{0k'}a_{k'})", \
     datafile."/vec_01k.dat" using 1:(abs($3)) with lines lw 2 lc 6 title "Im(Int p_{0k'}a_{k'})", \
     datafile."/vec_01k.dat" using 1:(abs($6)) with lines lw 2 lc 8 title "Re(Int p_{kk'}a_{k'})", \
     datafile."/vec_01k.dat" using 1:(abs($7)) with lines lw 2 lc 3 title "Im(Int p_{kk'}a_{k'})", \


#set output outfile2
#set ylabel 'b_{kw}(T)' font ', 15' enhanced
#unset yrange
#
#if (logscale == 1) {
#    set log y
#    plot datafile using 1:($3) with lines lw 2 lc 8 title '|b_{kw}(T)|^2'
#} else {
#    plot datafile using 1:3 with lines lw 2 title '|b_{kw}(T)|^2'
#}

unset log y

unset format y 
set output outfile2
plot datafile."/vec_01k.dat" using 1:($4)  pt 4 lc 4 title "Re(C^{/Symbol w}_{k0} p_{k0}a_0)", \
     datafile."/vec_01k.dat" using 1:($5)  pt 5 lc 6 title "Im(C^{/Symbol w}_{k0} p_{k0}a_0)", \
     datafile."/pk0.dat" using 1:($3) with lines lw 2 lc 8 title "p_{k0}", \
