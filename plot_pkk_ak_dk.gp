set terminal pngcairo size 1000,700
set output "pkl.png"

set title "Sum over k': pkk ak dk"
set xlabel "k"
set ylabel "value"

set grid
set key left top

plot \
"pkl.dat" using 1:2 with lines lw 2 title "pkk * ak dk", \
"pkl.dat" using 1:3 with lines dt 2 title "dkk * ak dk", \
"pkl.dat" using 1:4 with lines lw 3 title "analytic"


set terminal pngcairo size 1000,700
set output "pkl_imag.png"

set title "Imaginary part check"
set xlabel "k"
set ylabel "Im"

set grid

plot \
"pkl.dat" using 1:($2) with lines title "Im pkk", \
"pkl.dat" using 1:($3) with lines dt 2 title "Im dkk", \
"pkl.dat" using 1:($4) with lines lw 2 title "Im analytic"
