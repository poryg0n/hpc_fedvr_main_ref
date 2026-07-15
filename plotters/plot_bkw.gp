# Usage:
# gnuplot -c plot_bkw.gp datafile outfile

datafolder = ARG1
outfile  = ARG2
#mode     = int(ARG3)

#set format y "%.0s * 10^{%L}"
#unset format

set terminal pngcairo size 1000,700
set size ratio  0.3182 1,1


set datafile commentschars "#"


set output outfile."_bkw_bkwT.png"
set title "Momentum spectrum"
set xlabel "k (a.u.)"
set ylabel "b_{k{/Symbol w}}"
#set xrange [-1.6:1.6]


plot \
    datafolder."components_bkw.dat" using 2:($5) pt 6 dt (10,5) lw 2 lc 6 title "Re(b_{k{/Symbol w}})", \
    datafolder."components_bkw.dat" using 2:($6) pt 4 dt (10,5) lw 1 lc 2 title "Im(b_{k{/Symbol w}})", \


set output outfile."_hhg_bkw.png"
set title ""
set xlabel "k (a.u.)"
set ylabel "{/Symbol \362}|b_{k{/Symbol w}}|dk/(2pi)"
#set xrange [-1.6:1.6]

plot \
    datafolder."components_bkw2dk.dat" using 2:($3) with linesp pt 6 dt (10,5) lw 2 lc 6 title "|b_{k{/Symbol w}(T)}|^2", \
    datafolder."components_bkw2dk.dat" using 2:($5) with linesp pt 4 dt (10,5) lw 1 lc 2 title "|b_{k{/Symbol w}}|^2", \




set output outfile."_b0w_b0wT.png"
set title "Momentum spectrum"
set xlabel "{/Symbol w} (a.u.)"
set ylabel "b_{0{/Symbol w}}"
#set xrange [-1.6:1.6]

set grid
set key left 


plot \
    datafolder."components_b0w.dat" using 2:(($3)) with linesp pt 6 dt (10,2) lc 8 lw 2 title "Re(b_0{/Symbol w})", \
    datafolder."components_b0w.dat" using 2:(($4)) with linesp pt 6 dt (10,2) lc 6 lw 2 title "Im(b_0{/Symbol w})", \
    datafolder."components_b0w.dat" using 2:(($5)) with linesp pt 6 dt (10,9) lc 7 lw 2 title "Re(b_0{/Symbol w}(T))", \
    datafolder."components_b0w.dat" using 2:(($6)) with linesp pt 6 dt (10,9) lc 5 lw 2 title "Im(b_0{/Symbol w}(T))"




set format y "%2.0t{/Symbol \264} 10^{%L}"
set format y "10^{%L}"
set output outfile."_hhg_b0w.png"
set ylabel "|b_{0{/Symbol w}}|^2"
set xlabel "{/Symbol w} (a.u.)"
unset yrange
set xtics .5, .1, 1.1
set ytics 1.e-13, 1.e-2, 1.e2
set log y
plot \
    datafolder."components_b0w.dat" using 2:(($7**2)) with linesp pt 6 dt (10,5) lc 7 lw 2 title "|b_{0{/Symbol w}}|^2", \
    datafolder."components_b0w.dat" using 2:(($8**2)) with linesp pt 6 dt (10,5) lc 8 lw 2 title "|b_{0{/Symbol w}}(T)|^2", \
    datafolder.'/hhg.dat' using 1:2 with lines lw 2 title 'Eq.(104)'





