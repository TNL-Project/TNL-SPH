set terminal pngcairo enhanced dashed size 1000, 800 font "Verdana,25"
set output outputFileName

set xlabel "t(||{/:Bold g}||/H)^{1/2}"
set ylabel "p / (ρ||{/:Bold g}||H)"
#set title "Sensor P1"

set key top left
set key box maxcols 1

#set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
#set grid xtics lc rgb "#bbbbbb" lw 1 lt 0
set grid ytics lc rgb "black" lw 1 lt 0
set grid xtics lc rgb "black" lw 1 lt 0
set encoding utf8

#set xrange [2:10]
#set xrange [2:7]

# Normalisation
tCoef = ( 9.81 / 0.3 )**0.5 * 0.00002 * 100
pCoef = 1. / ( 1000. * 9.81 * 0.3 )

p num u ($1 * tCoef):($5 * pCoef) w l lc "black" lw 1 title "WCSPH-DBC", exp u 1:2 w l lt 2 lc "blue" lw 2  title " Lobovsky 2014",
