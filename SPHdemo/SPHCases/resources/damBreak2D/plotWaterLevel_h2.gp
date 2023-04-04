set terminal pngcairo enhanced dashed size 1000, 800 font "Verdana,30"
set output outputFileName

set xlabel "t(||{/:Bold g}||/H)^{1/2}"
set ylabel "h/H"
#set title "WaterLevel H2"

set key top left
set key box maxcols 1

#set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
#set grid xtics lc rgb "#bbbbbb" lw 1 lt 0
set grid ytics lc rgb "black" lw 1 lt 0
set grid xtics lc rgb "black" lw 1 lt 0
set encoding utf8

# Normalisation
tCoef = ( 9.81 / 0.3 )**0.5 * 0.00002 * 100
hCoef = 1. /  0.3

p num u ($1 * tCoef):($3 * hCoef) lc "black" ps 0.5 pt 5 title "WCSPH-DBC", exp u 1:2 lc "blue" ps 0.5 pt 5 title "Lobovsky 2014"
