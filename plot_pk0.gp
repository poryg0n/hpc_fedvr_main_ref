set terminal pngcairo size 1000,700
set output "pk0.png"

set title "pk0 consistency check"
set xlabel "k"
set ylabel "value"

set grid
set key left top

plot \
"pk0.dat" using 1:2 with lines lw 2 title "Re(pk0)", \
"pk0.dat" using 1:3 with lines dt 2 title "Re(p0k)", \
"pk0.dat" using 1:4 with lines dt 3 title "Re(dk0 -> p)", \
"pk0.dat" using 1:5 with lines lw 3 title "analytic"
