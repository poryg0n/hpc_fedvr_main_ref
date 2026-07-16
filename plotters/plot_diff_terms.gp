# Usage:
# gnuplot -c plot_pemd.gp datafolder output1 output2 logscale

datafolder = ARG1
outfile  = ARG2
logscale = int(ARG4)

set terminal pngcairo size 1000,700
set output outfile."_pkkak_1.png"

set xlabel 'k (a.u.)' font ', 15' enhanced
set ylabel '|C^{/Symbol w}_{k0} p_{k0}a_0|' font ', 15' enhanced
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
#    plot datafolder using 1:($2) with lines lw 2 lc 7 title 'C_{k{/Symbol w}}p_{k0}a_0'
#} else {
#    plot datafolder using 1:($2) with lines lw 2 lc 7 title 'C_{/Symbol w}_{k0}p_{k0}a_0'
#}

#
set ylabel '{/Symbol \362} C^{/Symbol w}_{k0} p_{0k}a_{k}dk' font ', 15' enhanced

 set log y
plot datafolder."/vec_01k.dat" using 1:(abs($2)) with linesp pt 4 dt (10,5) lw 2 lc 7 title "Re({/Symbol \362} C^{/Symbol w}_{k0} p_{0k}a_{k}dk)", \
     datafolder."/vec_01k.dat" using 1:(abs($3)) with linesp pt 5 dt (10,5) lw 2 lc 6 title "Im({/Symbol \362} C^{/Symbol w}_{k0} p_{0k}a_{k}dk)", \
     datafolder."/vec_01k.dat" using 1:(abs($6)) with linesp pt 6 dt (10,5) lw 2 lc 8 title "Re({/Symbol \362} C^{/Symbol w}_{kk'} p_{kk'}a_{k'}dk')", \
     datafolder."/vec_01k.dat" using 1:(abs($7)) with linesp pt 7 dt (10,5) lw 2 lc 3 title "Im({/Symbol \362} C^{/Symbol w}_{kk'} p_{kk'}a_{k'}dk')", \



set output outfile."_pkkak_2.png"

set ylabel 'C^{/Symbol w}_{k0} p_{k0}a_0' font ', 15' enhanced
set title 'Non integral term in b_{k{/Symbol w}}' font ",16" enhanced

unset log y
unset format
set grid
set size ratio 0.3182 1,1
set key left top
#set xrange [-1.33:1.33]
#set ytics  1.e-9, 1.e-2, 1.e-2
#set xtics -1.0, .5, 1.0
set yrange [-1.e-6 : 1.e-2]


plot datafolder."/vec_01k.dat" using 1:($2) with linesp pt 4 dt (10,5) lw 2 lc 7 title "Re({/Symbol \362} C^{/Symbol w}_{k0} p_{0k}a_{k}dk)", \
     datafolder."/vec_01k.dat" using 1:($3) with linesp pt 5 dt (10,5) lw 2 lc 6 title "Im({/Symbol \362} C^{/Symbol w}_{k0} p_{0k}a_{k}dk)", \
     datafolder."/vec_01k.dat" using 1:($6) with linesp pt 6 dt (10,5) lw 2 lc 8 title "Re({/Symbol \362} C^{/Symbol w}_{kk'} p_{kk'}a_{k'}dk')", \
     datafolder."/vec_01k.dat" using 1:($7) with linesp pt 7 dt (10,5) lw 2 lc 3 title "Im({/Symbol \362} C^{/Symbol w}_{kk'} p_{kk'}a_{k'}dk')", \




unset log y

unset format y 
unset yrange
set output outfile."_pk0a0.png"
plot datafolder."/vec_01k.dat" using 1:($4)  pt 4 lc 4 title "Re(C^{/Symbol w}_{k0} p_{k0}a_0)", \
     datafolder."/vec_01k.dat" using 1:($5)  pt 5 lc 6 title "Im(C^{/Symbol w}_{k0} p_{k0}a_0)", \
     datafolder."/pk0.dat" using 1:($3) with lines lw 2 lc 8 title "p_{k0}", \
