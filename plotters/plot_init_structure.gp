# Usage:
# gnuplot -c plot_pemd.gp datafile outputfile logscale

datafile = ARG1
outfile  = ARG2
logscale = int(ARG3)

if (!exists("file")) file = "fundamental.dat"



set terminal pngcairo size 1000,700
set output outfile


set title sprintf("Wavefunction: %s", file)
set xlabel "x"
set ylabel "Amplitude"
set grid

set format y "10^{%L}";
#set log y

plot \
    datafile using 2:4 w l lw 2 title "Re(ψ)", \
