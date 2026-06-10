set terminal pngcairo size 1000,700

datafile1 = ARG1
datafile2 = ARG2
datafile3 = ARG3
outfile1  = ARG4
outfile2  = ARG5
outfile3  = ARG6
label1    = ARG7
label2    = ARG8
label3    = ARG9

#mode     = int(ARG3)


set terminal pngcairo size 1000,700

set output outfile1
set ylabel "p_{kk'}a_{k'}{/Symbol D}k"

set size ratio  0.3182 1,1
#set xrange [.5:1.13]
#set yrange [1.e-14:1.e-2]
#set xtics .5, .1, 1.1
#set ytics 1.e-13, 1.e-2, 1.e-3

set format y "10^{%L}";

set title "|{/Symbol S}_{k'} p_{kk'} a_{k'}| dk"
set xlabel "k (a.u.)" font ", 15"
set ylabel "|{/Symbol S}_{k'} p_{kk'} a_{k'}| dk"  font ", 15"

set grid
set key center top
set key right bottom
set log y

set ylabel "|p_{kk'}a_{k'}|dk"
plot \
datafile1 using 1:2 with lines lw 2 title label1, \
datafile2 using 1:2 with lines lw 2 title label2, \
datafile3 using 1:2 with lines lw 3 title label3


set output outfile2
set ylabel "(E_k - E_{k'})|d_{kk'}a_{k'}|dk"
plot \
datafile1 using 1:4 with lines lw 2 title label1, \
datafile2 using 1:4 with lines lw 2 title label2, \
datafile3 using 1:4 with lines lw 3 title label3

set output outfile3
set ylabel "analytic"
plot \
datafile1 using 1:6 with lines lw 2 title label1, \
datafile2 using 1:6 with lines lw 2 title label2, \
datafile3 using 1:6 with lines lw 3 title label3


#set terminal pngcairo size 1000,700
#set output "pkl_imag.png"
#
#set title "Imaginary part check"
#set xlabel "k"
#set ylabel "Im"
#
#set grid
#
#plot \
#datafile using 1:($3) with lines title "Im pkk", \
#datafile using 1:($5) with lines dt 2 title "Im dkk", \
#datafile using 1:($7) with lines lw 2 title "Im analytic"
